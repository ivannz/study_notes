#! /usr/bin/env python
# -*- coding: UTF-8 -*-

__date__            = "2015-03-15"
__author__          = "Makhalova, Nazarov"
__email__           = "tpmakhalova@edu.hse.ru, innazarov@edu.hse.ru"
__dscription__      = """Основной вычислительный модуль"""

##########################################################################################
from matplotlib.figure import Figure
import matplotlib.pyplot as plotter
import numpy as np
import FillGaps as fg
import xlwt

from algorithm import initial as ini
from algorithm import criteria as crit
from algorithm import m_local as mloc

## Проверка динамического подключения модулей
def attach( module, procedure ) :
    try:
## ПОдключаем модуль
        return getattr( __import__( "algorithm." + module,
            globals( ), locals( ), [ procedure ], -1 ), procedure )
    except Exception, e:
## Пропускаем ошибку выше
        raise e

################################################################################
class model( object ):
    def __init__( self ):
        super( model, self ).__init__( )
        self.__dataset = self.__original = None
        self.__name = None
        self.__application = None
        self.__ready_for_clustering = False
        self.__current_partition = list( )
        self.__num_classes = 0
        self.__clust_criterion = self.__alpha = self.__p = self.__m_param = None
################################################################################
    def register( self, application ) :
        if self.__application is None :
            self.__application = application
        return self
    def start( self ) :
        pass
################################################################################
    def has_data( self ) :
        return self.__dataset is not None
    def load_datafile( self, fin ) :
        if fin is None : return
        self.__dataset = np.genfromtxt( fin.name, delimiter = ',' )
        self.__original = self.__dataset.copy( )
        self.__name = fin.name
        self.__ready_for_clustering = False
        self.__current_partition = list( )
        self.__actual_m = 1
        self.__xfer_count = 0
    def export_clustering_report( self, fp ) :
        if fp is None : return
# Object - class dictionary
        obj_class = {}
        for i in range(len(self.__current_partition)):
            cl = self.__current_partition[i]
            for j in range(len(cl)):
                obj_class[cl[j]]=i
## Write original data (from datafile)
        # with open('results.csv', 'w') as fp:
        n_rows = self.__original.shape[0]
        n_columns = self.__original.shape[1]
        # preparation to writing 
        font = xlwt.Font()
        font.bold = True
        alignment = xlwt.Alignment()
        alignment.horz = xlwt.Alignment.HORZ_CENTER
        header_style = xlwt.XFStyle() 
        header_style.font = font 
        header_style.alignment = alignment 
        header_style1 = xlwt.XFStyle() 
        header_style1.font = font 
        style_float = xlwt.XFStyle()
        style_float.num_format_str = '#,##0.000'
        
        # writeng to xls
        wb = xlwt.Workbook(encoding='utf-8')
        # Data output
        sheet0 = wb.add_sheet(u'Анализируемые данные')
        sheet0.write_merge(0, 0, 0, n_columns - 1, u'Исходные данные', header_style) 
        for i in range(n_rows):
            for j in range(n_columns):
                if np.isnan(self.__original[i][j]):
                    sheet0.write(i + 1, j, '')
                else:
                    sheet0.write(i + 1, j, self.__original[i][j],style_float)
        sheet0.write_merge(0, 0, n_columns + 1, 2 * n_columns, u'Данные без пропусков', header_style) 
        for i in range(n_rows):
            for j in range(n_columns):
                sheet0.write(i + 1, n_columns + 1 + j, self.__dataset[i][j],style_float)
        sheet0.write(0,2 * n_columns + 1, u'Номер кластера', header_style1)
        for i in range(n_rows):
            sheet0.write(i + 1 ,2 * n_columns + 1, obj_class[ i ] + 1)
        # Clustering summary
        sheet1 = wb.add_sheet(u'Результаты кластеризации')
        sheet1.write_merge(0, 0, 0, 10, u'Параметры кластеризации', header_style) 
        sheet1.write(1,0,u'Количество кластеров: {:1d}'.format(self.__num_classes))
        if (self.__clust_criterion == "sim_diff"):
            sheet1.write(2,0, u'Критерий качества разбиения: разность мер удаленности точек (внутри классов и вне)')
        else:
            sheet1.write(2,0, u'Критерий качества разбиения: отношение раззности мер удаленности точек (внутри классов и вне) к их сумме')      
        sheet1.write(3,0,u'Параметр альфа: {:.2f}'.format( self.__alpha))
        sheet1.write(4,0,u'Параметр p: {:.2f}'.format( self.__p))
        sheet1.write(5,0,u'Количество перебросов : {:d}'.format( self.__xfer_count ))
        sheet1.write(6,0,u'Максимальное количество точек при перебросе : {:d}'.format( self.__actual_m ))
        sheet1.write_merge(7, 8, 0, 0, u'Номер', header_style1)
        sheet1.write_merge(7, 8, 1, 1, u'Объём', header_style1)
        sheet1.write_merge(7, 7, 2, 2 + n_columns-1, u'Центроиды', header_style)
        sheet1.write_merge(7, 7, 2 + n_columns, 2 + 2*n_columns-1, u'ст. отклонение', header_style)
        for j in range( n_columns ) :
            sheet1.write(8, 2 + j, u'X%d'%j)
            sheet1.write(8, 2 + n_columns + j, u'X%d'%j)
        for i in range(len(self.__current_partition)):
            cl = self.__current_partition[ i ]
            cluster = self.__dataset[ cl ]
            sheet1.write(9 + i, 0, i+1)
            sheet1.write(9 + i, 1, len( cl ) )
            centroid = np.mean( cluster, axis = 0 )
            sdev = np.std( cluster-centroid, axis = 0 )
            for j in range( n_columns ) :
                sheet1.write(9 + i, 2 + j, centroid[ j ], style_float)
                sheet1.write(9 + i, n_columns + 2 + j, sdev[ j ], style_float)
        wb.save(fp.name) 

#		fp.writelines(u'Исходные данные\n')
#        for member in range(len(self.__original[0])):
#           fp.write('\tx{:d};'.format(member+1))
#        fp.write('\n')
#        for i in self.__original:
#            for j in i:
#               fp.write(str.format("{0:.3f}", j).replace('.',',')+';')
#            fp.write('\n')
#        fp.write('\n')
### Write results of clustering
#        fp.write('Параметры кластеризации\n'.encode('utf-8'))
#        fp.write(('Количество кластеров: {:1d}.\n'.format(self.__num_classes)).encode('utf-8'))
#        if (self.__clust_criterion == "sim_diff"):
#            fp.write(('Критерий качества разбиения: разность мер удаленности точек (внутри классов и вне).\n').encode('utf-8'))
#        else:
#            fp.write(('Критерий качества разбиения: отношение раззности мер удаленности точек (внутри классов и вне) к их сумме.\n').encode('utf-8'))
#        fp.write(('Параметр альфа: {:.2f}.\n'.format( self.__alpha)).encode('utf-8'))
#        fp.write(('Параметр p: {:.2f}.\n'.format( self.__p)).encode('utf-8'))
#        fp.write('Данные кластеризации (с заполненными пропусками)\n'.encode('utf-8'))
#        for i in range(len(self.__dataset)):
#            for j in self.__dataset[i]:
#                fp.write(str.format("{0:.3f}", j).replace('.',',')+';'+str(obj_class[i])+';')
#            fp.write('\n')
    def __call__( self ) :
        return self.__dataset
    def get_data_info( self ) :
        if self.__dataset is None : return ( '', 0, 0, )
        return ( self.__name, self.__dataset.shape[ 0 ], self.__dataset.shape[ 1 ] )
    def __prepare_data( self, original = False ) :
## Fill missing values if necessary
        data = self.__original.copy( )
        missing_mask = np.any( np.isnan( data ), axis = 1 )
        if not original :
            if np.any( missing_mask ) :
                data = fg.fill_missing( data )
## Prepare data for plottting
        if data.shape[ 1 ] > 3 :
            if original and np.any( missing_mask ) :
## Take the first three dimensions
                d3d = data[:,:3]
            else :
## Do the projection onto the 
                mu = np.mean( data, axis = 0 )
                sigma = np.std( data, axis = 0 )
                u, _, v = np.linalg.svd( ( data - mu ) / sigma )
                d3d = np.dot( data, v[:, :3] )
        else :
            d3d = np.zeros( ( data.shape[ 0 ], 3 ), np.float )
            d3d[ :, :data.shape[ 1 ] ] = data
        return d3d, missing_mask        
    def show_data( self, figure, original = False ) :
        if self.__dataset is None : return
## Plot the data
        d3d, missing_mask = self.__prepare_data( original )
        # ax = Axes3D( figure ) ;  ax.mouse_init( )
        ax = figure.gca( projection = '3d' )
        if len( self.__current_partition ) :
            # print self.__current_partition
            for cl, pi in zip( [ 'r', 'g', 'b', 'm', 'k', 'c', 'y' ], self.__current_partition ) :
                cur, c_m = d3d[ pi ], missing_mask[ pi ]
                ax.scatter( cur[np.logical_not( c_m ),0], cur[np.logical_not( c_m ),1], cur[np.logical_not( c_m ),2], marker = 'o', c = cl )
                ax.scatter( cur[c_m,0], cur[c_m,1], cur[c_m,2], marker = '^', c = cl )
        else :
## Highlight the datapoint which were reconsturcted
            kept_data = d3d[np.logical_not( missing_mask ) ]
            ax.scatter( kept_data[:,0], kept_data[:,1], kept_data[:,2], marker = 'o', c = 'b' )
            miss_data = d3d[ missing_mask ]
            ax.scatter( miss_data[:,0], miss_data[:,1], miss_data[:,2], marker = '^', c = 'r' )
        xLabel = ax.set_xlabel( '\nX', linespacing = 3.2 )
        yLabel = ax.set_ylabel( '\nY', linespacing = 3.1 )
        zLabel = ax.set_zlabel( '\nZ', linespacing = 3.4 )
        ax.dist = 10
## Commit to the device
    def get_cluster_info( self, num ) :
        if self.__dataset is None : return
        if len( self.__current_partition ) < 1 : return
        cur_inx = self.__current_partition[ num ]
        d3d, missing_mask = self.__prepare_data( False )
        cur, c_mask = d3d[ cur_inx ], missing_mask[ cur_inx ]
        return len( cur_inx ), len( c_mask[ c_mask ] ), np.mean( cur, axis = 0 )
    def show_cluster( self, figure, num ) :
        if self.__dataset is None : return
        if len( self.__current_partition ) < 1 : return
        d3d, missing_mask = self.__prepare_data( False )
        cur_inx = self.__current_partition[ num ]
        oth_inx = [ i for i in xrange( d3d.shape[ 0 ] ) if i not in cur_inx]
        cur, c_mask = d3d[ cur_inx ], missing_mask[ cur_inx ]
        oth, o_mask = d3d[ oth_inx ], missing_mask[ oth_inx ]
        ax = figure.gca( projection = '3d' )
## Plot other datapoints
        ax.scatter( oth[o_mask,0], oth[o_mask,1], oth[o_mask,2], marker = '^', c = '#D0D0D0' )
        ax.scatter( oth[np.logical_not(o_mask),0], oth[np.logical_not(o_mask),1], oth[np.logical_not(o_mask),2], marker = 'o', c = '#D0D0D0' )
## Plot the contents of this cluster
        ax.scatter( cur[c_mask,0], cur[c_mask,1], cur[c_mask,2], marker = '^', c = 'red' )
        ax.scatter( cur[np.logical_not(c_mask),0], cur[np.logical_not(c_mask),1], cur[np.logical_not(c_mask),2], marker = 'o', c = 'red' )
## Mark the centroid (reconsturcted data)
        ccur = np.mean( cur, axis = 0 )
        ax.scatter( ccur[0], ccur[1], ccur[2], marker = 'x', c = 'black', s = 50 )
        print cur_inx
        xLabel = ax.set_xlabel( 'X', linespacing = 3.2 )
        yLabel = ax.set_ylabel( 'Y', linespacing = 3.1 )
        zLabel = ax.set_zlabel( 'Z', linespacing = 3.4 )
        ax.dist = 10
## Commit to the device
        plotter.show( )
    def run_cluster( self ) :
        if not self.__ready_for_clustering :
            return
## Получаем данные из базы
        self.__dataset = self.__original.copy( )
## Пополняем пропуски
        self.__dataset = fg.fill_missing( self.__dataset )
## Выбираем начальное разбиение
        partition = ini.hierarchical( self.__dataset,
            num_classes = self.__num_classes, method = 'single' )
## Задаём функцию криетрия
        criterion = crit.criterion( self.__dataset, fun = self.__clust_criterion,
            alpha = self.__alpha, p = self.__p )
## Кластеризуем
        print "R0: ", partition
        print "I0: ", criterion.evaluate( partition )
        mlocal = mloc.m_local( criterion )
        self.__current_partition, status = mlocal( partition, m = self.__m_param )
        self.__actual_m, self.__xfer_count = status
        print "R1: ", partition
        print "I1: ", criterion.evaluate( partition )
        print "m,c: %d,%d"%(self.__actual_m, self.__xfer_count)
## Function to report on the lodad data and to set parameters
    def setup_begin( self ) :
        self.__ready_for_clustering = False
        pass
    def setup_end( self ) :
        self.__ready_for_clustering = True
        pass
    def get_avaliable_classes( self ) :
        if self.__dataset is None : return ( 0, 0, )
        classes = range( 2, self.__dataset.shape[ 0 ] // 15 + 1 )
        return classes, [ str( i ) for i in classes ]
    def get_avaliable_criteria( self ) :
        return [ "sim_diff", "sim_ratio", ], [ u"Разность", u"Отношение", ]
    def select_number_of_classes( self, choice ) :
        self.__num_classes = abs( int( choice ) )
    def read_alpha( self ) :
        return self.__alpha
    def read_p( self ) :
        return self.__p
    def read_number_of_classes( self ) :
        return self.__num_classes
    def select_criterion( self, choice ) :
        self.__clust_criterion = attach( "criteria", choice )
    def set_alpha( self, alpha ) :
        self.__alpha = abs( float( alpha ) )
    def set_p( self, p ) :
        self.__p = abs( float( p ) )
    def get_m_param_values( self ) :
        values = [ 1, 999 ]
        return values, [ str( v ) for v in values ]
    def set_m_param( self, m ) :
        self.__m_param = abs( int( m ) )
    def read_m_param( self ) :
        return self.__m_param

if __name__ == '__main__':
    pass