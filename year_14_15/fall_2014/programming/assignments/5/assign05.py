# Abstraction, Inheritance, Incapsulation, Polymoprhism
# Class, Object, Prototype

class person(object):
	def __init__(self, name):
		super(person, self).__init__()
		self.name = str(name)
	def __del__(self):
		pass
	def SayHello(self):
		print "Hi, my name is `%s`" % ( self.name )
	def trait(self):
		pass

class man(person):
	def __init__(self, name):
		super(man, self).__init__(name)
	def trait(self):
		print "I am a MAAAN!\n"

class woman(person):
	def __init__(self, name):
		super(woman, self).__init__(name)
	def trait(self):
		print "Tee-hee\n"


a = person("ivan")
a.SayHello()
a.trait()
b = man( 'ivan' )
b.trait()
c = woman( 'ivan' )
c.trait()

## Class categories
# Data layer
# Logic layer
# UI layer