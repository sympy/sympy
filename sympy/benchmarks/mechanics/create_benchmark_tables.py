import sqlite3

conn = sqlite3.connect('benchmarking.db')
c = conn.cursor()

# Table columns
table_cells = 'unix REAL, user TEXT, datestamp TEXT, platform TEXT, processor \
               TEXT, python_version TEXT, test_time REAL'

# Create the tables
c.execute('CREATE TABLE IF NOT EXISTS lagrange_1(%s)' % table_cells)
c.execute('CREATE TABLE IF NOT EXISTS kane_1(%s)' % table_cells)

# Close the cursor and database connection properly
c.close()
conn.close()
