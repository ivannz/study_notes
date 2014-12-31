#! /usr/bin/env python
# -*- coding: UTF-8 -*-

#### Python Programming Assignment 01 11/09/2014
## Написать консольную программу, которая на вход принимает начальную дикторию.
##   Она должна обходить файловую систему и сохранять полученную информацию в
##   текстовый файл, выделяя отступами папки и файлы. Программа также должна
##   замерять время обхода и сохранять его в заголовок файла.
## Выполнять первое задание и также хранить полученную информацию в графе.
##   Программа должна позволять добавлять и удалять папки и файлы из графа и
##   одновременно файловой системе, путем выделения текущей папки и текущего
##   файла.
## Программа должна позволять производить проверку на вхождение имени на
##   наличие в графе и выдавать информацию о том файл это или папка.

## Список комманд:
##   mount "путь к физической папке" "название точки"
##     -- загрузить фаловую структуру в память и прикрепить к указанной точке.
##   unmount "название точки"
##     -- Открепить присоединённую к указанной точке фаловую систему
##   ls [путь [шаблон]]
##     -- Отобразить содержимое текущей директории, или директории согласно
##     указанному относитеьному или абсолютному пути. При наличии параметра
##     "шаблон", являющегося регулярным выражением выводятся все папки и файлы,
##     чью имена удовлетворяют указанному щаблону.
##   cd [путь]
##     -- Отобразить текущую папку, если параметр "путь" не указан, или установить
##     текущкю папку в указанную относительным или абсолютным путём.
##   rm [y] [путь] шаблон
##     -- Удаляет все файлы и папки, чьи имена удовлетворяют шаблону, из папки,
##        указанной относительным или абсолютным путём. Первый необязательный
##        параметр определяет, следует ли требовать подтвержения для удаления.
##   rename "путь к объекту" "новое имя объекта"
##     --  Переименовывает указанный объект согласно новому имени.
##   save  "путь" "путь к физическому файлу" 
##     -- Сохраняет древовидную структуру объекта согласно указанному абсолютному
##        или относительному пути, по указанному физическому адресу в виде фала с
##        отступами.
##   set "флаг" on/off
##     -- Установить или определить указанный системный флаг в необходимое состояние.
##   get "флаг"
##     -- Отобразить состояние указанного системного флага. Доступные флаги:
##         virtual - установлен в состояние "ON" по умолчанию. Отвечает за псоледствия
##         от действий удаления или переименовывания на настоящую файловую систему,
##         загруженную в память. ИСПОЛЬЗОВАТЬ ОСТОРОЖНО!
##   exit
##     -- Выход
####################################################################################################
def ui_input(message) : return raw_input(message)

def ui_print(format, *args) : print format % args

def ui_error(format, *args) : print "ERROR: " + format % args

def ui_warning(format, *args) : print "Warning: " + format % args

def ui_prompt(message, choices = ['Y', 'N'], force = True) :
	options = set([opt.lower() for opt in list( choices )])
	prompt = "%s (%s): " % (message, "/".join(options))
	while True :
		answer = ui_input(prompt).lower()
		if answer.lower() in options or not force : break
		ui_print("Please make a valid choice!")
	return answer

####################################################################################################
import os
import re
import time
from assign01_node import inode as inode

def __is_type(node, kind) :
	return node.data.get("type", None) == kind.lower()

def __get_type(node) :
	return node.data.get("type", ".")

def __os_build_tree(path, name, parent) :
	if not name :
		return None
	path = os.path.expanduser(path)
	if not os.access(path, os.F_OK) :
		return None
	t_0 = time.clock()
	if os.path.isdir(path) :
		node = inode(name, parent, {"type": "d"})
		for name in os.listdir(path) :
			__os_build_tree(os.path.join(path, name), name, node)
	elif os.path.isfile(path) :
		node = inode(name, parent, {"type": "f"})
	else :
		return None
	stat = os.stat(path)
	node.data.update({"path": path, "size": stat.st_size, "scan": time.clock() - t_0})
	return node

def __os_remove(node) :
	if flag("virtual") :
		return True
	if __is_type(node, "f") :
		ui_print( "os.remove(\"%s\")", node.data["path"] )
	elif __is_type(node, "d") :
		ui_print( "os.rmdir(\"%s\")", node.data["path"] )
	return True

def __os_rename(node, new_name) :
	if flag("virtual") :
		return True
	if node.parent() is ROOT :
		return True
	if __is_type(node, "f") :
		ui_print( "os.renames(\"%s\", \"%s\")", node.data["path"], new_name )
	elif __is_type(node, "d") :
		ui_print( "os.rename(\"%s\", \"%s\")", node.data["path"], new_name )
	return True

####################################################################################################
CDIR = ROOT = inode("*", None, {"type":"d", "path":"*", "size":-1, "scan":-1.0})
FLAGS = { "Virtual": True, "Verbose": True }

def flag(flag, default=True) :
	return FLAGS.get(flag, default)

def set_cdir(node) :
	global CDIR
	CDIR = node

def __is_mounted(name) :
	return ROOT.fetch(name) is not None

def __unmount(name) :
	node = ROOT.fetch(name)
	if node is None :
		return None
	if node in CDIR.path() :
		set_cdir(ROOT)
	node.unlink()
	node.clear()
	return node

def __mount(name, real_path) :
	__unmount(name)
	return __os_build_tree(real_path, name, ROOT)

def __cd(node, path) :
	if node is None :
		node = ROOT
	if path is None :
		return node
	vertices = [vertex for vertex in path.split("/") if vertex and vertex != '.']
	if not vertices :
		return node
	if vertices[0] == ROOT.name() :
		node, vertices = ROOT, vertices[1:]
	for name in vertices :
		if not __is_type(node, "d") and node != ROOT :
			return node
		if name == ".." :
			if not node.is_root() :
				node = node.parent()
		else :
			next = node.fetch(name)
			if next is None :
				ui_print("Couldn't find \"%s\" in \"%s\"", name, node.uri())
				return None
			node = next
	return node

def __rename(node, new_name) :
	if node is not None and node != ROOT :
		old_name = node.rename(new_name)
		if old_name is not None :
			if __os_rename(node, new_name) :
				return node
			node.rename(old_name)
	return None

def __rm(node, pattern, prompt = True) :
	if node is None or node is ROOT:
		return None
	if node in CDIR.path() : set_cdir(node)
	children, undeleted = node.find( pattern ), set()
	for child in children :
		if prompt :
			reply = ui_prompt("Remove \"%s\"?" % child.uri(), ['y','n','all'])
			if reply == 'n' :
				undeleted.add(child)
				continue
			elif reply == 'all' : prompt = False
		## Line up the successors of the current node in a depth-first order,
		##  so that a parent is never visited before its children.
		vertices = set()
		for vertex, level in child.dfs() :
			if vertex in vertices :
				continue
			## What if one couldn't remove the real image?
			if not __os_remove( vertex ) :
				## Then the whole path from the current root to this
				##   unfortunate item must be kept!
				path = vertex.path( child )
				vertices.update(path[:-1])
		if not vertices :
			child.unlink()
			child.clear()
		undeleted.update(vertices)
	return undeleted

def __show(node) :
	if node is None :
		return None
	ui_print("Entity \"%s\"\nType: %s\nReal path:\"%s\"\nSize: %d\nScan time: %f",
		node.uri(), __get_type(node), node.data["path"], node.data["size"], node.data["scan"] )
	return node

def __save(node, real_path) :
	if node is None :
		return None
	if not real_path :
		real_path = "./assign01_test.txt"
	with open( real_path, "w" ) as fp :
		fp.write( "Header\n" )
		for child, level in node.dfs_tree() :
			fp.write(
				"\t" * level + child.name() + ( "/" if __is_type(child, "d") else "" ) +
				"\t" + str(child.data["scan"]) + ", " + "\"" + child.data["path"] + "\"" +
				"\n")
	return node

####################################################################################################ß
## Implement commands and functionality
## mount <real path> <drive>
def mount(args) :
	if len(args) != 2 :
		ui_print("Usage: mount <real path> <drive>")
		return False
	return __mount(args[1], args[0]) is not None

## unmount <drive>
def unmount(args) :
	if len(args) != 1 :
		ui_print("Usage: unmount <drive>")
		return False
	if not __is_mounted(args[0]) :
		ui_error("Drive \"%s\" was not found.", args[0])
		return False
	if flag("Verbose") :
		if ui_prompt("Un-mount virtual file system on \"%s\" (no actual files are affected)" % args[0] ) != 'y' :
			return True
	return __unmount(args[0]) is not None

## rm <pattern> [<path>] [<flags>]
def rm(args) :
	if len(args) < 1 :
		ui_print("Usage: rm [<flags>] [<path>] <pattern>")
		return False
	pattern = args[0]
	path = args[1] if len(args) > 1 else None
	flags = set([x.lower() for x in args[2:]])
	prompt = 'q' not in flags and flag("Verbose")
	undeleted = __rm(__cd(CDIR, path), pattern, prompt )
	if undeleted :
		for node in undeleted :
			ui_print("Couldn't remove \"%s\"", node.uri())
	return True

## rename <path> <new_name>
def ren(args) :
	if len(args) != 2 :
		ui_print("Usage: ren <path> <new_name>")
		return False
	node = __cd(CDIR, args[0])
	if node is None or node is ROOT :
		return False
	if flag("verbose") :
		if ui_prompt("Rename \"%s\" to \"%s\"?" % (node.uri(), args[1])) != 'y' :
			return False
	if __rename(node, args[1]) != node :
		ui_print( "Couldn't rename \"%s\" to \"%s\"", node, args[1] )
		return False
	return True

## info <path>
def info(args) :
	if len(args) != 1 :
		ui_print("Usage: info <path>")
		return False
	return __show(__cd(CDIR, args[0])) is not None

## save <path> <real destination>
def save(args) :
	if len(args) < 2 :
		ui_print("Usage: save <path> <real destination>")
		return False
	return __save(__cd(CDIR, args[0]), args[-1]) is not None

## ls [<path>] [<pattern>]
def ls(args) :
	path = args[0] if len(args) > 0 else None
	pattern = args[-1] if len(args) > 1 else None
	node = __cd(CDIR, path)
	if node is None :
		return False
	children = node.find( pattern )
	if len(children) < 1 :
		ui_print("Directory \"%s\" is empty", node.uri())
	else :
		if pattern :
			ui_print("Contents of \"%s\" matching \"%s\":", node.uri(), pattern)
		else :
			ui_print("Contents of \"%s\"", node.uri())
		for child in children :
			if __is_type(child, "d") :
				ui_print("\t-D\t%s", child.name())
			elif __is_type(child, "f") :
				ui_print("\t-F\t%s", child.name())
			else :
				ui_print("\t-.\t%s", child.name())
	return True

## cd [<path>]
def cd(args) :
	if len(args) > 0 :
		node = __cd(CDIR, args[0])
		if node :
			set_cdir( node )
	ui_print("Current directory is \"%s\".", CDIR)
	return True

## echo [<flags>] <message>
def echo(args):
	if len(args) < 1 :
		ui_print("Usage: echo <flags> <message>")
		return False
	flags, text = set([x.lower() for x in args[:-1]]), args[-1]
	for f in flags :
		if f == "u" :
			text = text.upper()
		elif f == "l" :
			text = text.lower()
	print(text)
	return True

## set <parameter> <on/off>
def set_flags(args) :
	global FLAGS
	if len(args) != 2 or args[1].lower() not in {"on", "off"} :
		ui_print("Usage: set <parameter> <on/off>")
		return False
	if flag("Verbose") and ui_prompt( "Are you sure?" ) != 'y' :
		return True
	FLAGS[args[0].lower()] = True if args[1].lower() == "on" else False
	return True

## get <parameter>
def get_flags(args) :
	if len(args) != 1 :
		ui_print("Usage: get <parameter>")
		return False
	if flag(args[0].lower(), None) is None :
		ui_print( "Parameter \"%s\" is undefined.", args[0] )
	else :
		ui_print( "Parameter \"%s\" is set to \"%s\"",
			args[0], "On" if FLAGS[args[0].lower()] else "Off")
	return True

## This dictionary stores the implemented functionality
commands = {
	"quit" 		: 	None,
	"mount" 	: 	mount,
	"unmount" 	: 	unmount,
	"cd" 		: 	cd,
	"ls" 		: 	ls,
	"rename" 	: 	ren,
	"rm" 		: 	rm,
	"set" 		: 	set_flags,
	"get" 		: 	get_flags,
	"show" 		: 	info,
	"save" 		: 	save,
	"echo" 		: 	echo
}

## The procedure, which parses the user input into a set of command line parameters
def parse_command(s):
	if not s : return (s,s,[])
## Split the input string by whitespace, unless it is within a double quoted string
	arg_list=[i for x in
		re.finditer(r"(?ux)\s*(?:\"((?:[^\\\"]|\\.)*)\"|(\S+))\s*", s)
			for i in x.groups() if i]
	return (s, arg_list[0].lower(), arg_list[1:])

## Main user interface loop
def ui() :
	ui_print("Python Disk Operating System\n -- Nazarov Ivan Spetember 2014")
	while True:
		cmd,command,args=parse_command(ui_input("%>"))
		if not command :
			continue
		elif command=="exit":
			break
		elif not commands.has_key(command):
			ui_error("Unknown command")
		else:
			fn = commands[command]
			if fn is None :
				ui_print("Command \"%s\" has no effect", cmd)
			elif not callable(fn) :
				ui_error("Internal error: object %s is badly implemented", cmd)
			else :
				if not fn(args) :
					ui_error("Errors during execution of %s", cmd)


################################################################################
ui()

# set virutal off
# y

# mount "./playground" "C:"

# ls

# cd "./C:"

# save "*/C:" "./playground/ASSIG01.TXT"

# rm .* "."

# exit

# """
# path = "*/Z:/R-CUDA"
# root = __cd(CDIR, path)

# os.chdir("/Users/user/Desktop/studies 2014-2015/programming/assign01/playground")

# os.mkdir("./a")
# open("./a", 'w').close()
# os.mkdir("./b")
# os.mkdir("./c")
# os.mkdir("./a/a")
# open("./a/a.txt", 'w').close()
# os.mkdir("./a/b")
# open("./a/b/a.txt", 'w').close()
# open("./a/b/b.txt", 'w').close()
# os.mkdir("./b/a")

# """
