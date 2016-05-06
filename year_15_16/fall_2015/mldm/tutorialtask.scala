// Linear sum
def linearSum( n : Long ) : Long = if( n > 0 ) n + linearSum( n - 1 ) else n

def linearSum_helper( n : Int, acc : Long ) : Long = {
    if( n  <= 0 ) {
        acc
    } else {
        linearSum2( n - 1, acc + n )
    }
}

def linearSum2( n : Int ) : Long = linearSum_helper( n, 0 )

def linearSum3( n : Int ) = Iterator.range( 1, n + 1 ) reduceLeft {_+_}

def linearSum4( n : Int ) = {
    var acc : Long = 0
    for { i <- Iterator.range( 1, n + 1 ) }{ acc = acc + i }
    acc
}

linearSum( 10 )
linearSum2( 10 )
linearSum3( 10 )
linearSum4( 10 )


def integrate( fn : Double => Double, interval : ( Double, Double ), n : Int = 10 ) = {
    val (a, b) = interval ; val delta : Double = ( b - a ) / n
    delta * ( Iterator.range( 0, n ) map {
            a + _ * delta
        } map {
            // x => fn( x + delta * 0.5 )
            x => 0.5 * ( fn( x ) + fn( x + delta ) )
        } reduceLeft {
            _ + _
        } )
}

integrate( math.sin( _ ), ( 0 , math.Pi ), 1000 )
