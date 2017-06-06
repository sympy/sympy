# This file will provide GUI access to display the benchmarking results

# Necessary module imports
import datetime as DT
import sqlite3
import tkinter as tk

# Import matplotlib backends
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
matplotlib.use('TkAgg')

# Establish a connection to the sqlite database
conn = sqlite3.connect('benchmarking.db')
c = conn.cursor()


# Create the tk app
class MainApp(tk.Tk):

    def __init__(self, *args, **kwargs):
        """This function will set up the tkinter GUI and will collect the list
        of tables present in the database"""

        # Initiate the base tk app code
        tk.Tk.__init__(self, *args, **kwargs)

        # Collect the names of the tables in the database
        c.execute("SELECT name FROM sqlite_master WHERE type='table';")
        self.tables = c.fetchall()

        # Set up a variable to keep track of the table currently being displayed
        self.tableviewing = tk.StringVar()
        self.tableviewing.set(self.tables[0][0])

        # Set up the tk menu that will allow switching between the different
        # tables in the database
        self.tablebutton = tk.Menubutton(self, text='Select Table',
                                         relief='raised')
        self.tablebutton.place(relx=0.2, rely=0.03, anchor='center')
        self.tablebutton.menu = tk.Menu(self.tablebutton)
        self.tablebutton['menu'] = self.tablebutton.menu
        for table in self.tables:
            self.tablebutton.menu.add_radiobutton(label=table[0],
                                                  variable=self.tableviewing,
                                                  value=table[0],
                                                  command=self.set_table)

        # Set up the button to switch between different python versions
        self.pyverbutton = tk.Menubutton(self, text='Select Python Version',
                                         relief='raised')
        self.pyverbutton.place(relx=0.5, rely=0.03, anchor='center')
        self.pyver = tk.StringVar()

        # Set up the button to switch between different platforms
        self.platbutton = tk.Menubutton(self, text='Select Platform',
                                        relief='raised')
        self.platbutton.place(relx=0.8, rely=0.03, anchor='center')
        self.plat = tk.StringVar()

        # Display the default table
        self.set_table()

    def filter_data(self):
        """This method takes the user selected filters and collects the desired
        data and graphs it"""

        # Reset the dates and times to empty lists
        self.dates = ['']
        self.times = []

        # Select the specific data values based on the selected filters
        if self.pyver.get() != 'All' and self.plat.get() != 'All':
            self.selected_data = [i for i in self.data if i[3] ==
                                  self.plat.get() and
                                  i[5].startswith(self.pyver.get())]
        elif self.pyver.get() != 'All':
            self.selected_data = [i for i in self.data if
                                  i[5].startswith(self.pyver.get())]
        elif self.plat.get() != 'All':
            self.selected_data = [i for i in self.data if i[3] ==
                                  self.plat.get()]
        else:
            self.selected_data = self.data

        # Create variables to graph dates and run times from the filtered data
        # list
        for row in self.selected_data:
            if DT.datetime.strptime(row[2], "%Y-%m-%d") == self.dates[-1]:
                self.times[-1] = (row[6] + self.times[-1])/2
            else:
                self.times.append(row[6])
                self.dates.append(DT.datetime.strptime(row[2], "%Y-%m-%d"))

        # Remove the blank entry from the dates list and graph the collected
        # rows
        self.dates.remove('')
        self.graph()

    def graph(self):
        """This method will graph the selected data and display it in the
        matplotlib canvas"""

        # Set up the matplotlib figure and provide the lables for its axes and
        # title
        frame = Figure(figsize=(5, 4), dpi=100)
        a = frame.add_subplot(111)
        a.plot(self.dates, self.times, self.dates, self.times, 'rs')
        a.set_xlabel('Dates')
        a.set_ylabel('Test Run Time (s)')
        a.set_title('Test Run Times for %s' %
                    self.tableviewing.get())

        # Put the matplotlib figure in the tkinter gui
        self.canvas = FigureCanvasTkAgg(frame, master=self)
        self.canvas.show()
        self.mpl_graph = self.canvas.get_tk_widget()
        self.mpl_graph.place(relx=0.5, rely=0.53, relheight=0.94, relwidth=1,
                             anchor='center')

    def readdata(self):
        """This method will collect all of the rows from the selected table
        along with the different platforms and python versions that the table
        contains"""

        # Collect all of the rows from the selected table
        c.execute('SELECT * FROM %s' % self.tableviewing.get())
        self.data = c.fetchall()

        # Reinitialize the python versions and platforms lists and add the
        # default all arguement
        self.pyvers = set(['All'])
        self.plats = set(['All'])

        # Collect the python versions and platforms available in this specific
        # table
        for row in self.data:
            self.pyvers.add(row[5][0:3])
            self.plats.add(row[3])

    def set_table(self):
        """This method resets the filters and the different python version and
        platform options when switching to view a different table"""

        # Reset the python version and platform filters when switching tables
        self.pyver.set('All')
        self.plat.set('All')

        # Collect and display the data
        self.readdata()
        self.filter_data()

        # Add the python versions and platforms to their respective menu buttons
        self.pyverbutton.menu = tk.Menu(self.pyverbutton)
        self.pyverbutton['menu'] = self.pyverbutton.menu
        for ver in self.pyvers:
            self.pyverbutton.menu.add_radiobutton(label=ver,
                                                  variable=self.pyver,
                                                  value=ver,
                                                  command=self.filter_data)

        self.platbutton.menu = tk.Menu(self.platbutton)
        self.platbutton['menu'] = self.platbutton.menu
        for ver in self.plats:
            self.platbutton.menu.add_radiobutton(label=ver,
                                                 variable=self.plat,
                                                 value=ver,
                                                 command=self.filter_data)


# Run the tk app
if __name__ == "__main__":
    app = MainApp()
    app.geometry('800x600')
    app.title("Benchmarking Viewer")
    app.mainloop()

    # Close the cursor and database connection properly
    c.close()
    conn.close()
