import datetime as dt

import sqlite3 as lite
con = lite.connect( "/Volumes/elements/avito_sqlite/database.sqlite" )
cur = con.cursor(  )

import pandas as pd
res = cur.execute( """ select * from searchinfo where userid = 444291 ; """).fetchall( )

user_searches = pd.DataFrame( res, columns = [ "SearchID", "SearchDate", "UserID", "IsUserLoggedOn", "IPID", "SearchQuery", "SearchLocationID", "SearchCategoryID", "SearchParams", ] )
user_searches[ "SearchDate" ] = user_searches[ "SearchDate" ].map( lambda x: dt.datetime.strptime( x, "%Y-%m-%d %H:%M:%S.0" ) )


import os

PATH = os.path.realpath( os.path.join( "study_notes", "data_study" ) )
user_searches.sort( "SearchDate" ).to_csv( os.path.join( PATH, "444291.csv" ), encoding = "utf8", index = False )
