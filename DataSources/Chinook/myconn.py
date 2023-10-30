# http://dev.mysql.com/doc/connector-python/en/index.html
#
import mysql.connector as mc
from mysql.connector import errorcode 
try:
  cnx = mc.connect(user='root', database='test')
  cursor = cnx.cursor()

  query = ("SELECT uname, fname, email, phone FROM info")
  cursor.execute(query)

  for (uname, fname, email, mobile) in cursor:
    print("%s : %s : %s : %s\n" % (uname,fname,email,mobile) )

  cursor.close()
  cnx.close

except mc.Error as err:
  print(err)

