
/* Variables and constants */
val a1 = 2        // constant
var b1 : Int = 1  // variable
b1 += 1

/* functions */
def add( x : Int , y : Int ) : Int = {
  x + y
}

def inc( x : Int ) = x + 1

/* lambda anonymous functions */
{ ( x:Int, y:Int ) =>  val z = x+y ; z }
{ ( x:Int, y:Int ) =>  x + y  }

/* functions : scope and closure */
var a : Int = 10;
val f = { x:Int => x+a }
f( 9 )
a = 1
f( 1 )

def func( a : Int ) = {
  { x : Int => x + a }
}

func(1)(2)

def func( b : Int ) : Int => Int = {
  { x : Int => x + b + a }
}

a = 1
func(10)(20)
a = 9
func(10)(20)
// I see you like functions, so I put a function inside a function, that returns a function

def epsilon( f : Int => Int, x : Int )  = f( x )

a = 100
epsilon( f, 10 )

a = 10
epsilon( f, 10 )

/* "def" declares a method, whereas "val" declares a field */

/* collections*/

// List are single-linked lists. cdr and car ( first, rest )
List(1,2,3)
4 :: List(1,2,3)
List(1,2,3) :: List(1,2,3)

val xs = List(1,2,3,4,5)
// takes a function and applies it to the whole collection.
xs.map( { x: Int => x + 10 } )
xs.map( { _ + 100 } )

xs.map( { x: Int => List( x, x*2, x*3 ) } )

val ys = xs.flatMap( { x: Int => List( x, x*2, x*3 ) } )

// conditional: evaluates the condition and then evaluate the relevant branch
if (ys %% 2 == 0) { throw "1" } else { x + 1 }

ys.filter( { x => x%2 == 0 } )

// Reducers: foldLeft
ys.foldLeft( 0 ) { _ + _ }
ys.foldLeft( 0 ) { (x: Int, y: Int) => x + y}

// Reducers: foldLeft
ys.foldRight( 0 ) { (x: Int, y: Int) => x + y}

ys.foldRight( "" ) { _ + _.toString() }
ys.foldLeft( "" ) { _ + _.toString() }

// Reducers: reduce (can't change types)
ys.reduceLeft { _ - _ }
ys.reduceRight { _ - _ }

// collections : sets
val xl = Set( 1,1,3,4,5 )

xl.map( { _ + 1 } )
xl.map { _ + 1 }
xl map { _  / 2 }

xl filter { _ != 3 }

xl filter { x => if(x % 2 == 0) true else false }

// Iteration
for { x <- ys ; if x %2 == 0 } yield x + 1
ys filter( _ % 2 == 0 ) map { _ + 1 }

for { x <- ys ; y <- xl ; if x % 2 == 0 } yield x + 1 - y

// Recursion
def g( n : Int, x : Int ) : Int = if ( n > 1 ) g( n - 1, x ) + x else x

g(10, 2)

def g( n : Int ) : Int = if ( n > 1 ) g( n - 1 ) + n else n
g(10)

//
val x = (1,3)
val a, b = x
val (u, v) = x


// A standalone scirpt file
object main extends Nothing {
	def def main(args: Array[String]): Unit = {
	  println( args )
	}
	
}

