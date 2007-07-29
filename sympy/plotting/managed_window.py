from pyglet.gl import *
from pyglet.window import Window
from pyglet.clock import Clock

from threading import Thread, Lock, Event

gl_lock = Lock()

class ManagedWindow(Window):

    fps_limit = 40
    default_win_args = dict(width=400,
                            height=300,
                            vsync=False,
                            resizable=True)

    def __init__(self, **win_args):
        self.win_args = dict(self.default_win_args, **win_args)
        Thread(target=self.__event_loop__).start()

    def __event_loop__(self, **win_args):
        gl_lock.acquire()
        try:
            try:
                super(ManagedWindow, self).__init__(**self.win_args)
                self.switch_to()
                self.setup()
            except Exception, e:
                print "Window Setup Error: %s" % str(e)
        finally:
            gl_lock.release()

        clock = Clock()
        clock.set_fps_limit(self.fps_limit)
        while not self.has_exit:
            dt = clock.tick()
            gl_lock.acquire()
            try:
                try:
                    self.switch_to()
                    self.dispatch_events()
                    self.clear()
                    self.update(dt)
                    self.draw()
                    self.flip()
                except Exception, e:
                    print "Window Event Loop Error: %s" % str(e)
            finally:
                gl_lock.release()
        super(ManagedWindow, self).close()

    def close(self):
        self.has_exit = True

    def setup(self):
        pass

    def update(self, dt):
        pass

    def draw(self):
        pass
