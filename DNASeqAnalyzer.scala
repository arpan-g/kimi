import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.SparkConf
import org.apache.spark.Partitioner
import sys.process._
import scala.sys.process._

import java.io._
import java.nio.file.{Paths, Files}
import java.net._
import java.util.Calendar
import java.util.Arrays
import scala.sys.process.Process
import scala.io.Source
import scala.collection.JavaConversions._
import scala.util.Sorting._

import net.sf.samtools._

import tudelft.utils.ChromosomeRange
import tudelft.utils.DictParser
import tudelft.utils.Configuration
import tudelft.utils.SAMRecordIterator

import org.apache.commons.lang3.exception.ExceptionUtils
import org.apache.spark.storage.StorageLevel._
import org.apache.spark.HashPartitioner

import collection.mutable.HashMap

object DNASeqAnalyzer 
{



final val MemString = "-Xmx14336m" 
final val RefFileName = "ucsc.hg19.fasta"
final val SnpFileName = "dbsnp_138.hg19.vcf"
final val ExomeFileName = "gcat_set_025.bed"


//////////////////////////////////////////////////////////////////////////////

def getListOfFiles(dir: String):List[File] = {
        val d = new File(dir)
        if (d.exists && d.isDirectory) {
          d.listFiles.filter(_.isFile).toList
        } else {
          List[File]()
        }
      }


def bwaRun (files:List[File],sc:SparkContext,config: Configuration) :org.apache.spark.rdd.RDD[(Int, net.sf.samtools.SAMRecord)] =
{
	//var kvPairs = Array[(Int, SAMRecord)]()
	var sa = 0
	var fileList = List[File]()
	for(file <- files if (!file.getName.contains(".crc")&& (!file.getName.contains("_SUCCESS")))){
	if(file.exists){
		fileList ::= file
		}
	}
	val fileRDD = sc.parallelize(fileList,config.getNumInstances.toInt)
	val exitValue = fileRDD.mapPartitions({iter => for(i<-iter) yield createSamFile(i)},true)
	exitValue.take(20).foreach(println)
	val resFiles = getListOfFiles("/data/result")
	val resRDD = sc.parallelize(resFiles,config.getNumInstances.toInt)
	val kvPairs= resRDD.mapPartitions({iter => for(i<-iter) yield createBWAKVPairs(i)},true)
	val kvFlat = kvPairs.flatMap(x=>x)
        return kvFlat
}
def writeToBAM(fileName: String, samRecordsSorted: Array[SAMRecord], config: Configuration) : ChromosomeRange = 
{
	val header = new SAMFileHeader()
	header.setSequenceDictionary(config.getDict())
	val outHeader = header.clone()
	outHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
	val factory = new SAMFileWriterFactory();
	val writer = factory.makeBAMWriter(outHeader, true, new File(fileName));
	
	val r = new ChromosomeRange()
	val input = new SAMRecordIterator(samRecordsSorted, header, r)
	while(input.hasNext()) 
	{
		val sam = input.next()
		writer.addAlignment(sam);
	}
	writer.close();
	
	return r
}
def createBWAKVPairs(file:File):Array[(Int, SAMRecord)] =
{
	val bwaKeyValues = new BWAKeyValues(file.getAbsolutePath)
	bwaKeyValues.parseSam()
	val kvPairs: Array[(Int, SAMRecord)] = bwaKeyValues.getKeyValuePairs()
	return kvPairs
}

def createSamFile(file:File):String=
{
val cmd = Seq("/data/spark/tools/bwa","mem","-p", "-t","8","/data/spark/ref/ucsc.hg19.fasta",file.getAbsolutePath)
	cmd #> new File("/data/result/"+file.getName+".sam") !
	
return file.getName
}

def variantCall (chrRegion: Int, samRecordsSorted: Array[SAMRecord], config: Configuration ) : Int =
{	
	val tmpFolder = config.getTmpFolder
	val toolsFolder = config.getToolsFolder
	val refFolder = config.getRefFolder
	val numOfThreads = config.getNumThreads
	Arrays.sort(samRecordsSorted,new SAMCompare())
	
	// Following is shown how each tool is called. Replace the X in regionX with the chromosome region number (chrRegion). 
	// 	You would have to create the command strings (for running jar files) and then execute them using the Scala's process package. More 
	// 	help about Scala's process package can be found at http://www.scala-lang.org/api/current/index.html#scala.sys.process.package.
	//	Note that MemString here is -Xmx14336m, and already defined as a constant variable above, and so are reference files' names.
	
	// SAM records should be sorted by this point
	val path = tmpFolder+"/region"+chrRegion+"-p1.bam"
	val chrRange = writeToBAM(path, samRecordsSorted, config)
	// Picard preprocessing
	val pic1= Seq("java",MemString,"-jar",toolsFolder+"/CleanSam.jar" ,"INPUT="+path,"OUTPUT="+tmpFolder+"/region"+chrRegion+"-p2.bam") !
		
	val pic2 = Seq("java",MemString,"-jar",toolsFolder+"/MarkDuplicates.jar","INPUT="+tmpFolder+"/region"+chrRegion+"-p2.bam","OUTPUT="+tmpFolder+"/region"+chrRegion+"-p3.bam", 
			"METRICS_FILE="+tmpFolder+"/region"+chrRegion+"-p3-metrics.txt") !
	val pic3 = Seq("java", MemString,"-jar",toolsFolder+"/AddOrReplaceReadGroups.jar","INPUT="+tmpFolder+"/region"+chrRegion+"-p3.bam","OUTPUT="+tmpFolder+"/region"+chrRegion+".bam", 
			"RGID=GROUP1", "RGLB=LIB1", "RGPL=ILLUMINA", "RGPU=UNIT1", "RGSM=SAMPLE1") !
	val pic4 = Seq("java",MemString,"-jar",toolsFolder+"/BuildBamIndex.jar","INPUT="+tmpFolder+"/region"+chrRegion+".bam") !
	
	//val rm1 = Seq("rm",tmpFolder+"/region"+chrRegion+"-p1.bam",tmpFolder+"/region"+chrRegion+"-p2.bam",tmpFolder+"/region"+chrRegion+"-p3.bam",tmpFolder+"/region"+chrRegion+"-p3-metrics.txt") !
	
	// Make region file 
		val tmpBed = new File(tmpFolder+"/tmp"+chrRegion+".bed")
		chrRange.writeToBedRegionFile(tmpBed.getAbsolutePath())
		val pic5 = Seq(toolsFolder+"/bedtools", "intersect", "-a", refFolder+"/"+ExomeFileName, "-b",tmpFolder+"/tmp"+chrRegion+".bed","-header")
		pic5 #> new File(tmpFolder+"/bed"+chrRegion+".bed") !
	//	val rm2 = Seq("rm",tmpFolder+"/tmp"+chrRegion+".bed") !
	
	// Indel Realignment 
		val pic6 =Seq("java", MemString, "-jar", toolsFolder+"/GenomeAnalysisTK.jar","-T","RealignerTargetCreator", "-nt",numOfThreads, "-R", refFolder+RefFileName,"-I", tmpFolder+"/"+"region"+chrRegion+".bam","-o",tmpFolder+"/region"+chrRegion+".intervals","-L",tmpFolder+"/bed"+chrRegion+".bed") !
		val pic7 = Seq("java",MemString, "-jar", toolsFolder+"/GenomeAnalysisTK.jar","-T","IndelRealigner","-R",refFolder+RefFileName, "-I", tmpFolder+"/region"+chrRegion+".bam","-targetIntervals",tmpFolder+"/region"+chrRegion+".intervals","-o",tmpFolder+"/region"+chrRegion+"-2.bam","-L",tmpFolder+"/bed"+chrRegion+".bed") !
	//	delete tmpFolder/regionX.bam, tmpFolder/regionX.bai, tmpFolder/regionX.intervals
//val rm3 = Seq("rm",tmpFolder+"/region"+chrRegion+".bam",tmpFolder+"/region"+chrRegion+".bai",tmpFolder+"/region"+chrRegion+".intervals") !
// Base quality recalibration 
val pic8 = Seq("java",MemString, "-jar",toolsFolder+"/GenomeAnalysisTK.jar","-T","BaseRecalibrator","-nct", numOfThreads, "-R",refFolder+RefFileName,"-I",tmpFolder+"/region"+chrRegion+"-2.bam","-o", tmpFolder+"/region"+chrRegion+".table","-L",tmpFolder+"/bed"+chrRegion+".bed","--disable_auto_index_creation_and_locking_when_reading_rods","-knownSites", refFolder+SnpFileName) !
val pic9 = Seq("java",MemString,"-jar",toolsFolder+"/GenomeAnalysisTK.jar","-T","PrintReads","-R", refFolder+RefFileName,"-I",tmpFolder+"/region"+chrRegion+"-2.bam","-o", tmpFolder+"/region"+chrRegion+"-3.bam","-BQSR",tmpFolder+"/region"+chrRegion+".table", "-L",tmpFolder+"/bed"+chrRegion+".bed") !
//	val rm4 = Seq("rm",tmpFolder+"/region"+chrRegion+"-2.bam", tmpFolder+"/region"+chrRegion+"-2.bai", tmpFolder+"/region"+chrRegion+".table") !
	
	// Haplotype -> Uses the region bed file
val pic10 = Seq("java", MemString,"-jar",toolsFolder+"/GenomeAnalysisTK.jar", "-T", "HaplotypeCaller", "-nct", numOfThreads,"-R", refFolder+RefFileName ,"-I",tmpFolder+"/region"+chrRegion+"-3.bam", "-o", tmpFolder+"/region"+chrRegion+".vcf", "-stand_call_conf", "30.0", "-stand_emit_conf","30.0","-L", tmpFolder+"/bed"+chrRegion+".bed","--no_cmdline_in_header","--disable_auto_index_creation_and_locking_when_reading_rods") !
	
	// delete tmpFolder/regionX-3.bam, tmpFolder/regionX-3.bai, tmpFolder/bedX.bed
	//val rm4 = Seq(rm,tmpFolder+"/region"+chrRegion+"-3.bam", tmpFolder+"/region"+chrRegion+"-3.bai", tmpFolder+"/bed"+chrRegion+".bed") !
	// return the content of the vcf file produced by the haplotype caller.
	//	Return those in the form of <Chromsome number, <Chromosome Position, line>>*/

	return pic10
}


def interleave(input:Iterator[(Long, (List[String], List[String]))]): Iterator[String] =
{
 var res = List[String]()
while(input.hasNext){	
  val next = input.next
  val aiter = next._2._1.toIterator
  val biter = next._2._2.toIterator
 
  while (aiter.hasNext && biter.hasNext)
  { 
    res ::= aiter.next
 if(aiter.hasNext)	
    res ::= aiter.next
 if(aiter.hasNext)
    res ::= aiter.next
 if(aiter.hasNext)
    res ::= aiter.next
 if(biter.hasNext)
    res ::= biter.next
 if(biter.hasNext)
    res ::= biter.next 
 if(biter.hasNext)
    res ::= biter.next 
 if(biter.hasNext)
    res ::= biter.next 
  }
}
  res.toIterator
}

def main(args: Array[String]) 
{
	val config = new Configuration()
	config.initialize()
		 
	val conf = new SparkConf().setAppName("DNASeqAnalyzer")
	// For local mode, include the following two lines
	conf.setMaster("local[" + config.getNumInstances() + "]")
	conf.set("spark.cores.max", config.getNumInstances())
	
	// For cluster mode, include the following commented line
	//conf.set("spark.shuffle.blockTransferService", "nio") 
	
	val sc = new SparkContext(conf)
	val fast1q = sc.textFile("/data/spark/fastq/fastq1.fq")
	val fast2q = sc.textFile("/data/spark/fastq/fastq2.fq")
	val fast1q_with_index_reduceByKey = fast1q.zipWithIndex().map{case(s,index)=>(index/4,List(s))}
	val fast1q_with_index = fast1q_with_index_reduceByKey.reduceByKey((a,b)=>(b:::a))
	val fast2q_with_index_reduceByKey = fast2q.zipWithIndex().map{case(s,index)=>(index/4,List(s))}
	val fast2q_with_index = fast2q_with_index_reduceByKey.reduceByKey((a,b)=>(b:::a))
	val joined_map = fast1q_with_index.join(fast2q_with_index,new ExactPartitioner(config.getNumInstances.toInt))
	val joinpartition = joined_map.mapPartitions(iter=>interleave(iter))
	println(joinpartition.count)
	joinpartition.saveAsTextFile(config.getInputFolder)
	val files = getListOfFiles(config.getInputFolder)
	val kvPairs = bwaRun (files,sc,config)
	val kvChrRegion = kvPairs.map(a=>((a._1-1)/3,a._2))
	val kvGroup =	kvChrRegion.groupByKey(new RegionPartition())
	val kvPartition = kvGroup.mapPartitions({iter => for(i<-iter) yield variantCall(i._1,i._2.toArray,config)},true)
	kvPartition.count	
	val ifiles = getListOfFiles(config.getTmpFolder())

	var fileList = List[File]()
	for(file <- ifiles if (file.getName.contains(".vcf")&& (!file.getName.contains("idx")))){
		if(file.exists){
			fileList ::= file
		}
	}
	
	var fileRDD: org.apache.spark.rdd.RDD[String]  = null
	for(file <- fileList ){	
		if(file.exists){
		val ran = sc.textFile(file.getAbsolutePath)
		if(fileRDD!=null)
		fileRDD = fileRDD.union(ran)
		else
		fileRDD = ran	
		}	
	}
	val filterOutTemp = fileRDD.filter(!_.contains("##"))
	filterOutTemp.count
	val filterOut = filterOutTemp.filter(!_.contains("#"))
	filterOut.count
	val splitMap=filterOut.map(a=>(a.split("\t"),a))
	val abc = splitMap.map(a=>(a._1(0).replaceAll("chr","").replaceAll("X","23").toInt->(a._1(1).toInt,a._2)))	 
	val a = abc.sortByKey(true)
	a.count
	a.take(20).foreach(println)
	val pw = new PrintWriter(new File("/home/agovindaraju/Documents/big-data/spark-1.5.0-bin-hadoop2.4/samRecord.vcf"))
	pw.write("#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ls\n")
	a.collect().foreach(ln=>pw.write(ln._2._2+"\n"))
	pw.close()
// Rest of the code goes here
}
//////////////////////////////////////////////////////////////////////////////
} // End of Class definition

// /data/spark/tools/bwa mem -p -t 8 /data/spark/ref/ucsc.hg19.fasta






