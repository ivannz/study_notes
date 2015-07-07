#! /usr/bin/env python
# -*- coding: UTF-8 -*-

## Common dependencies
import numpy as __np

class m_local(object) :
	"""Класс m-локальной оптимизации разбиения. Использование:
		sloc = m_local( crit = объект критерия качества кластеризации )
	Методы: полный раунд s-локальной оптимизации. Разбиение изменяется
	или не изменяетя автоматически. Возвращает новое значение критерия
	кластеризации для полученного нового разбиения. Агрумент partition
	передаётся по ссылке и претерпевает изменения!
		sloc(
			partition = текущее разбиение,
			s = количество перебрасываемых точек )
	Полная процедура m-локальной оптимизации. Производит m полных
	раундов s-локальной оптимизации. Аргумент разбиения передаётся
	по ссылке и изменяется в процессе работы процедуры.
		sloc(
			partition = текущее разбиение,
			m = количество перебрасываемых точек )
	"""
	def __init__( self, criterion ):
		super(m_local, self).__init__( )
		self.__criterion = criterion
## A basic step of the s-local optimisation procedure
	def __xfer( self, partition, src, dst, V0, s = 1, _xfer_count = 0 ) :
## There is no point in shuffling the data points withing the same class
		if src == dst :
			return V0, _xfer_count
## Keep moving batches of s elements partition until convergence
		while len( partition[ src ] ) >= s + 2 :
## Find candidates for moving
			sNN = self.__criterion.find_candidates( partition[ src ], partition[ dst ], s = s )
## Save the original classes (copies) 
			S0, D0 = list( partition[ src ] ), list( partition[ dst ] )
## Move them from S to D: this modifies the partition partition directly!
			list.extend( partition[ dst ], ( S0[ n ] for n in sNN ) )
			partition[ src ] = list( x for n, x in enumerate( S0 ) if n not in sNN )
			_xfer_count += 1
## Compute the criterion of the modified partition
			V1 = self.__criterion.evaluate( partition )
## If the criterion has not increased, rollback
			if V1 <= V0 :
				partition[ src ], partition[ dst ] = S0, D0
				break
## Indicate that the partition has been altered and proceed
			V0 = V1
## The partition is updated (or not) automatically
		return V0, _xfer_count
	def __cycle( self, partition, src, V0, s = 1, _xfer_count = 0 ) :
## a full cycle through the classes
		dst = 0
## Pick the class to modify. If we fall off the array, this means stabilization has occurred
		while dst < len( partition ) and len( partition[ src ] ) >= s + 2 :			
			V1, _xfer_count = self.__xfer( partition, src, dst, V0,
				s = s, _xfer_count = _xfer_count )
## Continue if the transfers did not yield significant fittness increase
			if V1 <= V0 :
				dst += 1
				continue
## Otherwise restart
			V0 = V1 ; dst = 0
		return V0, _xfer_count
	def s_local( self, partition, V0 = None, s = 1, _xfer_count = 0 ) :
		if V0 is None :
			V0 = self.__criterion.evaluate( partition )
		if len( partition ) < 2 :
			return V0, _xfer_count
		src = 0
		while src < len( partition ) :
## If the class is of sufficient volume ( the least volume is 2 points )
			if len( partition[ src ] ) >= s + 2 :
##  begin the s-local point transfer cycle
				V1, _xfer_count = self.__cycle( partition, src,
					V0, s = s, _xfer_count = _xfer_count )
## If the partition has been modified, restart
				if V1 > V0 :
					V0 = V1 ; src = 0
					continue
## Otherwise continue
			src += 1
## The s-local optimization step is finished when no
##  modification to the partition have been made.
		return V0, _xfer_count
## m-local optimization: go through the number of points transferred between
##  classes in succession.
	def __call__( self, partition, m = 5 ) :
		V0 = self.__criterion.evaluate( partition )
		if len( partition ) < 2 :
			return V0, _xfer_count
		s = 1 ; max_s = 1 ; _xfer_count = 0
		while s <= m :
			V1, _xfer_count = self.s_local( partition,
				V0, s = s, _xfer_count = _xfer_count )
			if V1 <= V0 :
## No improvement -- proceed to the next round of point transfers
				s += 1
				continue
## Remember the last size of a bundle before restarting
			if max_s < s :
				max_s = s
## Restart from the very beginning
			V0 = V1 ; s = 1
		return partition, (max_s, _xfer_count)
