import java.io.File
import java.io.InputStream
import java.io.FileInputStream
import java.util._
import net.sf.samtools.util.BufferedLineReader
import net.sf.samtools._
import tudelft.utils._

class BWAKeyValues(filePath: String)
{
	val keyValues = scala.collection.mutable.ArrayBuffer.empty[(Int, SAMRecord)]
	val mFile = new File(filePath);
	val is = new FileInputStream(mFile);
    val validationStringency: SAMFileReader.ValidationStringency = SAMFileReader.ValidationStringency.LENIENT;
    val mReader = new BufferedLineReader(is);
    val samRecordFactory = new DefaultSAMRecordFactory();
	private var mCurrentLine: String = null
	
	def getKeyValuePairs() : Array[(Int, SAMRecord)] = 
	{
		return keyValues.toArray
	}
	
    def writePairedSAMRecord(sam: SAMRecord) : Integer = 
	{
        var count = 0
        val read1Ref = sam.getReferenceIndex()
        val read2Ref = sam.getMateReferenceIndex()
		
		if (!sam.getReadUnmappedFlag() && (read1Ref > 0) && (read1Ref <= 24))
		{
			var region = read1Ref
		
			keyValues.append((region, sam))
			count = count + 1;
			
		}
		
		return count
    }
		
	def advanceLine() : String = 
    {
        mCurrentLine = mReader.readLine()
        return mCurrentLine;
    }
	
	def parseSam() =  
	{
		var mParentReader: SAMFileReader = null
        val headerCodec = new SAMTextHeaderCodec();
        headerCodec.setValidationStringency(validationStringency)
        val mFileHeader = headerCodec.decode(mReader, mFile.toString())
        val parser = new SAMLineParser(samRecordFactory, validationStringency, mFileHeader, mParentReader, mFile)
        // now process each read...
        var count = 0
        mCurrentLine = mReader.readLine()
		
        while (mCurrentLine != null) 
		{
			val samrecord = parser.parseLine(mCurrentLine, mReader.getLineNumber())
			count += writePairedSAMRecord(samrecord)
			//advance line even if bad line
			advanceLine();
		}
        
       // System.out.println("SAMstream counts " + count + " records");
    }
}
