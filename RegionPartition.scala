import org.apache.spark.Partitioner

class RegionPartition[V]()  extends Partitioner {

   def numPartitions = 8
  def getPartition(key: Any): Int = {
    val k = key.asInstanceOf[Int]
    return k
  }
}
