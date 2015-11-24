import org.apache.spark.Partitioner

class ExactPartitioner[V](partitions: Int)  extends Partitioner {

  def numPartitions = partitions
  def getPartition(key: Any): Int = {
    val k = key.asInstanceOf[Long]
    val par = (k.toInt % partitions)
    // `k` is assumed to go continuously from 0 to elements-1.
    return par
  }
}
