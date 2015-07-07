#! /usr/bin/env python
# -*- coding: UTF-8 -*-

__date__            = "2015-03-16"
__author__          = "Makhalova, Nazarov"
__email__           = "tpmakhalova@edu.hse.ru, innazarov@edu.hse.ru"
__status__          = "Alpha"
__version__         = "0.9"
__dscription__      = """Основной модуль работы по курсу "Структурно-классификационные методы интеллектуального анализа данных и прогнозирования в слабо формализованных системах" """

import sys
reload(sys)
sys.setdefaultencoding( 'utf-8' )
import time as tm

try :
    import tkinter as tk
    from tkinter.ttk import Combobox, Entry, Button
    from tkinter import filedialog as file_dlg
except :
    import Tkinter as tk
    from ttk import Combobox, Entry, Button
    import tkFileDialog as file_dlg

##########################################################################################
class Combox( Combobox ):
    def __init__( self, values, labels, **kwargs ) :
        self.values, self.labels = values, labels
        Combobox.__init__( self, values = labels, **kwargs )
    def current( self ) :
        return self.values[ Combobox.current( self ) ]

##########################################################################################
class Application( tk.Frame ):
    def __init__( self, hWnd, model ) :
        # super( Application, self ).__init__( )
        tk.Frame.__init__( self, hWnd )
        self.option_add( '*tearOff', False )
        self.__wnd = hWnd
        self.__menubar = None
        self.__menuitems = {}
        self.__model = model.register( self )
    def start( self ) :
## Parent
        self.__wnd.title( "Анализ данных" )
        self.__wnd.geometry( '{}x{}'.format( 600, 150 ) )
        self.__wnd.resizable( False, False )
## Menu bar
        self.__menubar = tk.Menu( self.__wnd, tearoff = 0 )
        self.__wnd.config( menu = self.__menubar )
## Exit
        # self.__menubar.add_command( label = "Выход", underline = 0, command = self.__cmd_menu_exit )
## Data
        self.__menuitems[ 'data' ] = m_data = tk.Menu( self.__menubar )
        self.__menubar.add_cascade( label = "Данные", underline = 0, menu = m_data )
        m_data.add_command( label = "Загрузить данные", command = self.__cmd_menu_data_open )
        m_data.add_command( label = "Сохранить данные", command = self.__cmd_menu_data_save )
        m_data.entryconfig( 1, state = tk.DISABLED )
## Show
        self.__menuitems[ 'show' ] = m_show = tk.Menu( self.__menubar )
        self.__menubar.add_cascade( label = "Просмотр", underline = 0, menu = m_show )
        # self.__menubar.entryconfig( 2, state = tk.DISABLED )
        m_show.add_command( label = "Сырые данные", command = self.__cmd_menu_show_view_original )
        m_show.add_command( label = "Данные без пропусков", command = self.__cmd_menu_show_view_processed )
        m_show.add_separator( )
        m_show.add_command( label = "Результаты", command = self.__cmd_menu_show_results )
        # m_show.entryconfig( 0, state = tk.DISABLED )
        # m_show.entryconfig( 1, state = tk.DISABLED )
        # m_show.entryconfig( 3, state = tk.DISABLED )
## clustering : Add to the main bar
        self.__menuitems[ 'clustering' ] = m_clust = tk.Menu( self.__menubar )
        self.__menubar.add_cascade( label = "Кластеризация", underline = 0, menu = m_clust )
        # self.__menubar.entryconfig( 3, state = tk.DISABLED )
        m_clust.add_command( label = "Запуск", command = self.__cmd_menu_cluster_run )
        # m_clust.entryconfig( 0, state = tk.DISABLED )
## About : Add to the main bar
        self.__menuitems[ 'about' ] = m_about = tk.Menu( self.__menubar )
        self.__menubar.add_cascade( label = "Помощь", underline = 0, menu = m_about )
        m_about.add_command( label = "о программе", command = self.__show_about_dialog )
        m_about.add_command( label = "руководство пользователя", command = self.__show_manual )
## Initialize the controller
        self.__model.start( )
## Invoke the dispatcher
        self.__wnd.mainloop( )
## View -- window command routines
    def __cmd_menu_exit( self ) :
        self.quit( )
    def __cmd_menu_data_open( self ) :
        fin = self.__show_open_dialog( )
        self.__model.load_datafile( fin )
        self.__display_show_datafile( )
        fin.close( )
        self.__menuitems[ 'data' ].entryconfig( 1, state = tk.ACTIVE )
    def __cmd_menu_data_save( self ) :
        fin = self.__show_save_dialog( )
        self.__model.export_clustering_report( fin )
    def __cmd_menu_show_view_original( self ) :
        plt = figure_window( tk.Toplevel( self ), title = u"Исходные данные", modal = True )
        fig = plt.figure( figsize = ( 8, 6 ), facecolor = 'w', dpi = 90 )
        self.__model.show_data( fig, original = True )
    def __cmd_menu_show_view_processed( self ) :
        plt = figure_window( tk.Toplevel( self ), title = u"Данные без пропусков", modal = True )
        fig = plt.figure( figsize = ( 8, 6 ), facecolor = 'w', dpi = 90 )
        self.__model.show_data( fig, original = False )
    def __cmd_menu_show_results( self ) :
        result_window( tk.Toplevel( self ), self.__model ).start( )
    def __cmd_menu_cluster_run( self ) :
        clustering_window( tk.Toplevel( self ), self.__model ).start( )
    def __display_show_datafile( self ) :
        if not self.__model.has_data( ) : return
## Show basic info on the loaded datafile
        filename, n, attr = self.__model.get_data_info( )
        tk.Label( self.__wnd, text = u"Загружены данные из файла %s" % filename ).grid( row = 0, sticky = tk.W )
        tk.Label( self.__wnd, text = u"Количество объектов: %d" % n ).grid( row = 1, sticky = tk.W )
        tk.Label( self.__wnd, text = u"Количество признаков: %d" % attr ).grid( row = 2, sticky = tk.W )
## Enable menu options
        # self.__menuitems[ 'show' ].entryconfig( 0, state = tk.ACTIVE )
        # self.__menuitems[ 'show' ].entryconfig( 1, state = tk.ACTIVE )
        # self.__menubar.entryconfig( 2, state = tk.ACTIVE )
        # self.__menubar.entryconfig( 3, state = tk.ACTIVE )
    def __display_error( self, error ) :
        err_wnd = tk.Toplevel( self )
        err_wnd.geometry( '{}x{}'.format( 300, 40 ) )
        err_wnd.resizable( False, False )
        tk.Label( err_wnd, text = error ).grid( row = 0, sticky = tk.W )
    def __show_open_dialog( self ) :
        return file_dlg.askopenfile(
            filetypes = ( ( "CSV", "*.csv" ), ( "All files", "*.*" ) ) )
    def __show_save_dialog( self ) :
        return file_dlg.asksaveasfile( mode = 'w',
			filetypes = ( ( u"Книга Excel 97-2003", "*.xls" ), ( "All files", "*.*" ) ) )
    def __show_about_dialog( self ) :
        about = about_window( tk.Toplevel( self ), title = u"о программе", modal = True )
        about.show( )
    def __show_manual( self ) :
        import os
        import win32com.client as win32
        if os.name == 'nt' :
            word = win32.gencache.EnsureDispatch( 'Word.Application' )
            word.Documents.Open( os.path.join( os.path.realpath( '.' ),
                u"Руководство пользователя программы анализа данных.docx" ) )
            word.Visible = True


##########################################################################################
##########################################################################################
class about_window( tk.Frame ):
    def __init__(self, hWnd, title, modal = False ):
        tk.Frame.__init__( self, hWnd )
        hWnd.title( title )
        hWnd.geometry( '{}x{}+50+50'.format( 363, 120 ) )
        hWnd.resizable( False, False )
        if modal:
            hWnd.grab_set( )
        hWnd.bind( '<Escape>', self.__close )
        self.__wnd = hWnd
    def show( self, **kwargs ) :
## The number of classes
        tk.Label( self.__wnd, text = "Авторы:" ).grid( row = 1, column = 0, sticky = tk.W )
        tk.Label( self.__wnd, text = __email__ ).grid( row = 1, column = 1, sticky = tk.W )
        T = tk.Text( self.__wnd, height = 9, width = 45 )
        T.grid( row = 2, column = 0, columnspan = 2, sticky = tk.W )
        T.insert( tk.END, __dscription__ + u"\nАвторы: " + __email__ )
    def __close( self, event ) :
        self.__wnd.destroy( )

##########################################################################################
##########################################################################################
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
class figure_window( tk.Frame ):
    def __init__(self, hWnd, title, modal = False ):
        tk.Frame.__init__( self, hWnd )
        hWnd.title( title )
        hWnd.resizable( False, False )
        if modal:
            hWnd.grab_set( )
        hWnd.bind( '<Escape>', self.__close )
        self.__wnd = hWnd
        self.__figure = None
    def figure( self, **kwargs ) :
        self.__figure = Figure( **kwargs )
        canvas = FigureCanvasTkAgg( self.__figure, master = self.__wnd )
        canvas._tkcanvas.pack( side = tk.TOP, fill = tk.BOTH, expand = 1 )
        canvas.show( )
        ax = Axes3D( self.__figure )
        ax.mouse_init( )
        return self.__figure
    def __close( self, event ) :
        self.__wnd.destroy( )

##########################################################################################
##########################################################################################
class result_window( tk.Frame ):
    def __init__(self, hWnd, model, callback = None ):
        tk.Frame.__init__( self, hWnd )
        self.__wnd = hWnd
        self.__callback = None
        self.__model = model
        self.start( )
    def start( self ) :
        if not self.__model.has_data( ) :
            self.__wnd.destroy( )
            return
        self.__model.setup_begin()
        self.__wnd.title( u"Результаты кластеризации" )
        self.__wnd.geometry( '{}x{}+50+50'.format( 380, 120 ) )
        self.__wnd.resizable( False, False )
        self.__wnd.grab_set( )
        self.__wnd.bind( '<Escape>', self.__close )
## The number of classes
        cls_num = self.__model.read_number_of_classes( )
        classes = range( cls_num ) ; labels = [ str( c + 1 ) for c in classes ]
        self.cur_class = Combox( classes, labels, master = self.__wnd, width = 5, state = 'readonly' )
        tk.Label( self.__wnd, text = "Класс:" ).grid( row = 1, column = 0, sticky = tk.W )
        self.cur_class.grid( row = 1, column = 1, sticky = tk.W )
        self.cur_class.set( labels[-1] )
        self.cur_class.bind( '<<ComboboxSelected>>', self.__onChoice )
    def __close( self, event ) :
        self.__wnd.destroy( )
    def __onChoice( self, event ) :
        title = u"Класс #{}".format( self.cur_class.current( ) + 1 )
        size, miss, cent = self.__model.get_cluster_info( self.cur_class.current( ) )
        tk.Label( self.__wnd, text = u"Точек в классе %d" % size ).grid( row = 2, column = 0, sticky = tk.W )
        tk.Label( self.__wnd, text = u"Реконструировано точек %d" % miss ).grid( row = 3, column = 0, sticky = tk.W )
        plt = figure_window( tk.Toplevel( self ), title = title, modal = False )
        fig = plt.figure( figsize = ( 8, 6 ), facecolor = 'w', dpi = 90 )
        self.__model.show_cluster( fig, self.cur_class.current( ) )

##########################################################################################
##########################################################################################
class clustering_window( tk.Frame ):
    def __init__(self, hWnd, model, callback = None ):
        tk.Frame.__init__( self, hWnd )
        self.__wnd = hWnd
        self.__callback = None
        self.__model = model
        self.start( )
    def start( self ) :
        if not self.__model.has_data( ) :
            self.__wnd.destroy( )
            return
        self.__model.setup_begin()
        self.__wnd.title( u"Кластеризации" )
        self.__wnd.geometry( '{}x{}+50+50'.format( 380, 150 ) )
        self.__wnd.resizable( False, False )
        self.__wnd.bind( '<Escape>', lambda e: self.onClose( ) )
        self.__wnd.protocol( "WM_DELETE_WINDOW", self.onClose  )
        self.__wnd.grab_set( )
## The number of classes
        cls_id, cls_label = self.__model.get_avaliable_classes( )
        self.num_classes = Combox( cls_id, cls_label , master = self.__wnd, width = 5, state = 'readonly' )
        tk.Label( self.__wnd, text = "Количество классов:"
            ).grid( row = 1, column = 0, sticky = tk.W )
        self.num_classes.grid( row = 1, column = 1, sticky = tk.W )
        cls_num = self.__model.read_number_of_classes( )-min(cls_id)
        if cls_num < 1 :
            cls_num = -1
        self.num_classes.set( cls_label[ cls_num ] )
## The target criterion
        crit_id, crit_label = self.__model.get_avaliable_criteria( )
        self.crit_fun = Combox( crit_id, crit_label, master = self.__wnd, width = 12, state = 'readonly' )
        tk.Label( self.__wnd, text = "Критерий качества:"
            ).grid( row = 3, column = 0, sticky = tk.W )
        self.crit_fun.grid( row = 3, column = 1, sticky = tk.W )
        self.crit_fun.set( crit_label[-1] )
## The similarity matrix parameters
        self.alpha_box = Entry( master = self.__wnd, width = 5 )
        tk.Label( self.__wnd, text = "Параметр Альфа:"
            ).grid( row = 4, column = 0, sticky = tk.W )
        self.alpha_box.grid( row = 4, column = 1, sticky = tk.W )
        alpha = self.__model.read_alpha( )
        self.alpha_box.insert( 0, ".5" if alpha is None else str( alpha ) )
        self.p_box = Entry( master = self.__wnd, width = 5 )
        tk.Label( self.__wnd, text = "Параметр P:"
            ).grid( row = 5, column = 0, sticky = tk.W )
        self.p_box.grid( row = 5, column = 1, sticky = tk.W )
        p = self.__model.read_p( )
        self.p_box.insert( 0, "8" if p is None else str( p ) )
## The optimisation parameter
        m_par_id, m_par_label = self.__model.get_m_param_values( )
        self.m_param = Combox( m_par_id, m_par_label , master = self.__wnd, width = 5, state = 'normal' )
        tk.Label( self.__wnd, text = "Параметр m-локальной оптимизации:"
            ).grid( row = 6, column = 0, sticky = tk.W )
        self.m_param.grid( row = 6, column = 1, sticky = tk.W )
        m_par = self.__model.read_m_param( )
        self.m_param.set( m_par_label[-1] if m_par is None else str( m_par ) )
## Add a button to start
        self.submit = Button( master = self.__wnd, width = 12,
            state = tk.ACTIVE, text = u"Запуск", command = self.onClose )
        self.submit.grid( row = 7, column = 1, sticky = tk.W )
    def onClose( self ) :
        self.__model.select_number_of_classes( self.num_classes.current( ) )
        self.__model.select_criterion( self.crit_fun.current( ) )
        self.__model.set_alpha( self.alpha_box.get( ) )
        self.__model.set_p( self.p_box.get( ) )
        self.__model.set_m_param( self.m_param.current( ) )
        self.__model.setup_end( )
        self.__model.run_cluster( )
        self.__wnd.destroy( )

##########################################################################################
if __name__ == '__main__' :
    from model import model as mdl
    Application( tk.Tk( ), mdl( ) ).start( )
    sys.exit( 0 )
