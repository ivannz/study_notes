package io.crayfis.apps
import org.apache.spark._
import org.apache.spark.rdd._

// val sconfig = new SparkConf( )
// sconfig.setMaster( "local[*]" )

// val sc = new SparkContext( sconfig )

val lines = sc.textFile( "file:///Users/user/4300-8.txt" )

// val words = lines flatMap { _.split( " " ) }
val words = lines flatMap { line : String => line.split( " " ) }
// val total = words map { w:String => 1 } sum

val freqs = words map { ( _, 1 ) } reduceByKey { _ + _ }
// freqs take 10 mkString "\n"

// val total = words map { ( _, 1 ) } map { w, f => f } sum
def wordDictionary( ws: Iterator[ String ] ) : Iterator[ ( String, Int ) ] = {
    ws map { (_,1) } // reduceByKey { _ + _ }
}

val ns = words.mapPartitions( wordDictionary ).reduceByKey { _ + _ } map { case (w, freq) => freq }

val total = ns sum( )

val freq = ns map{ n => n.toDouble / total }

.foldLeft( "" ) { _ + _.toString() }