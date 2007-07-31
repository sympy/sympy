from pyglet.window import key
from pyglet.window.mouse import LEFT, RIGHT
from util import get_direction_vectors, get_basis_vectors

class PlotController(object):
    
    normal_mouse_sensitivity = 4.0
    modified_mouse_sensitivity = 1.0

    normal_key_sensitivity = 160.0
    modified_key_sensitivity = 40.0

    keymap = {
                key.LEFT:'left',
                key.A:'left',
                key.NUM_4:'left',
                
                key.RIGHT:'right',
                key.D:'right',
                key.NUM_6:'right',
                
                key.UP:'up',
                key.W:'up',
                key.NUM_8:'up',
                
                key.DOWN:'down',
                key.S:'down',
                key.NUM_2:'down',
                
                key.Z:'rotate_z_neg',
                key.NUM_1:'rotate_z_neg',
                
                key.C:'rotate_z_pos',
                key.NUM_3:'rotate_z_pos',
                
                key.Q:'spin_left',
                key.NUM_7:'spin_left',
                key.E:'spin_right',
                key.NUM_9:'spin_right',

                key.X:'reset_rotations',
                key.NUM_5:'reset_rotations',
                
                key.NUM_ADD:'zoom_in',
                key.PAGEUP:'zoom_in',
                key.R:'zoom_in',

                key.NUM_SUBTRACT:'zoom_out',
                key.PAGEDOWN:'zoom_out',
                key.F:'zoom_out',
                
                key.RSHIFT:'modify_sensitivity',
                key.LSHIFT:'modify_sensitivity',
             }

    def __init__(self, window):
        self.window = window
        self.action = {
                # Rotation around the view Y (up) vector
                'left':False,
                'right':False,
                # Rotation around the view X vector
                'up':False,
                'down':False,
                # Rotation around the view Z vector
                'spin_left':False,
                'spin_right':False,
                # Rotation around the model Z vector
                'rotate_z_neg':False,
                'rotate_z_pos':False,
                # Reset to the default rotation
                'reset_rotations':False,
                # Performs camera z-translation
                'zoom_in':False, 
                'zoom_out':False,
                # Use alternative sensitivity (speed)
                'modify_sensitivity':False,
            }
        
    def update(self, dt):
        z = 0
        if self.action['zoom_out']: z -= 1
        if self.action['zoom_in']: z += 1
        if z != 0:
            self.window.camera.zoom_relative(z/10.0, self.get_key_sensitivity()/10.0)        
        
        dx, dy, dz = 0, 0, 0
        if self.action['left']: dx -= 1
        if self.action['right']: dx += 1
        if self.action['up']: dy -= 1
        if self.action['down']: dy += 1
        if self.action['spin_left']: dz += 1
        if self.action['spin_right']: dz -= 1
        if dx != 0:
            self.window.camera.euler_rotate(dx*dt*self.get_key_sensitivity(), *(get_direction_vectors()[1]))
        if dy != 0:
            self.window.camera.euler_rotate(dy*dt*self.get_key_sensitivity(), *(get_direction_vectors()[0]))
        if dz != 0:
            self.window.camera.euler_rotate(dz*dt*self.get_key_sensitivity(), *(get_direction_vectors()[2]))

        rz = 0
        if self.action['rotate_z_neg']: rz -= 1
        if self.action['rotate_z_pos']: rz += 1
        if rz != 0:
            self.window.camera.euler_rotate(rz*dt*self.get_key_sensitivity(), *(get_basis_vectors()[2]))

        if self.action['reset_rotations']:
            self.window.camera.init_rot_matrix()

        return True

    def get_mouse_sensitivity(self):
        if self.action['modify_sensitivity']:
            return self.modified_mouse_sensitivity
        else:
            return self.normal_mouse_sensitivity

    def get_key_sensitivity(self):
        if self.action['modify_sensitivity']:
            return self.modified_key_sensitivity
        else:
            return self.normal_key_sensitivity

    def on_key_press(self, symbol, modifiers):
        if symbol in self.keymap:
            self.action[self.keymap[symbol]] = True

    def on_key_release(self, symbol, modifiers):
        if symbol in self.keymap:
            self.action[self.keymap[symbol]] = False

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        if buttons & LEFT:
            self.window.camera.spherical_rotate((x-dx,y-dy),(x,y), self.get_mouse_sensitivity())
        if buttons & RIGHT:
            self.window.camera.translate(dx, dy, self.get_mouse_sensitivity())

    def on_mouse_scroll(self, x, y, dx, dy):
        self.window.camera.zoom_relative(dy, self.get_mouse_sensitivity())
