import kivy
kivy.require('1.5.1')
from kivy.uix.button import Button
from kivy.uix.gridlayout import GridLayout
from kivy.uix.textinput import TextInput
from sympy import *
from kivy.app import App
from kivy.uix.label import Label
from kivy.core.window import Window

class entrywidget(GridLayout):

    global flag
    global l
    l = Label()
    flag = 0


    def _keyboard_close(self):
        print('My keyboard have been closed!')
        self._keyboard = None

    def __init__(self, **kwargs):
        self._keyboard = Window.request_keyboard(self._keyboard_close, self)
        if self._keyboard.widget:
            vkeyboard = self._keyboard.widget
            vkeyboard.layout = 'numeric'
        super(entrywidget,self).__init__(**kwargs)
        self.cols = 2
        self.add_widget(Label(text = 'Expression'))
        self.text_input = TextInput(multiline = False)
        self.add_widget(self.text_input)
        self.add_widget(Label(text = "Result"))
        self.text_input.bind(on_text_validate = self.on_enter)

    def on_enter(self, value):
        global flag
        global l
        print flag
        if(flag == 1):
            self.remove_widget(l)
        z = self.text_input.text
        arr = z.split('=')
        x = Symbol("x")
        y = Symbol("y")
        z = Symbol("z")
        t = Symbol("t")
        if(len(arr) > 1):
            x = eval(arr[1])
            t = eval(arr[1])
            t = str(t)
        else:
            t = eval(z)
            t = str(t)
            l = Label(text = t)
            self.add_widget(l)
            flag = 1

		
class SympyApp(App):
    def build(self):
        return entrywidget()

if __name__ == '__main__':
    SympyApp().run()

