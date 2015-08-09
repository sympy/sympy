from sympy import Symbol
from sympy.physics.vector import Point, ReferenceFrame
from sympy.physics.mechanics import inertia, Body


# To test if two frames were created the same way.
def check_reference_frames(frame1, frame2):
    assert frame1.name == frame2.name
    assert frame1.str_vecs == frame2.str_vecs
    assert frame1.indices == frame2.indices
    assert frame1.varlist == frame2.varlist
    assert frame1._dlist == frame2._dlist
    assert frame1._ang_vel_dict == frame2._ang_vel_dict
    assert frame1._ang_acc_dict == frame2._ang_acc_dict
    assert frame1._dcm_dict == frame2._dcm_dict


def test_default():
    body = Body('body')
    assert body.name == 'body'
    assert body.force_list == []
    point = Point('body_masscenter')
    point.set_vel(body.frame, 0)
    com = body.masscenter
    frame = body.frame
    assert com.vel(frame) == point.vel(frame)
    assert body.mass == Symbol('body_mass')
    assert body.inertia == (inertia(body.frame, 1, 1, 1), body.masscenter)


def test_custom_rigid_body():
    # Body with RigidBody.
    rigidbody_masscenter = Point('rigidbody_masscenter')
    rigidbody_mass = Symbol('rigidbody_mass')
    rigidbody_frame = ReferenceFrame('rigidbody_frame')
    body_inertia = inertia(rigidbody_frame, 1, 0, 0)
    rigid_body = Body('rigidbody_body', rigidbody_masscenter, rigidbody_mass,
               rigidbody_frame, body_inertia)
    com = rigid_body.masscenter
    frame = rigid_body.frame
    rigidbody_masscenter.set_vel(rigidbody_frame, 0)
    assert com.vel(frame) == rigidbody_masscenter.vel(frame)
    assert com.pos_from(com) == rigidbody_masscenter.pos_from(com)

    assert rigid_body.mass == rigidbody_mass
    assert rigid_body.inertia == (body_inertia, rigidbody_masscenter)


def test_particle_body():
    #  Body with Particle
    particle_masscenter = Point('particle_masscenter')
    particle_mass = Symbol('particle_mass')
    particle_frame = ReferenceFrame('particle_frame')
    particle_body = Body('particle_body', particle_masscenter, particle_mass,
                  particle_frame)
    com = particle_body.masscenter
    frame = particle_body.frame
    particle_masscenter.set_vel(particle_frame, 0)
    assert com.vel(frame) == particle_masscenter.vel(frame)
    assert com.pos_from(com) == particle_masscenter.pos_from(com)

    assert particle_body.mass == particle_mass
    assert not hasattr(particle_body, "_inertia")


def test_particle_body_add_force():
    #  Body with Particle
    particle_masscenter = Point('particle_masscenter')
    particle_mass = Symbol('particle_mass')
    particle_frame = ReferenceFrame('particle_frame')
    particle_body = Body('particle_body', particle_masscenter, particle_mass,
                  particle_frame)

    a = Symbol('a')
    particle_body.add_force((a, 0, 0), particle_body.masscenter)
    assert len(particle_body.force_list) == 1
    point = particle_body.masscenter.locatenew(
        particle_body._name + '_point0', 0)
    point.set_vel(particle_body.frame, 0)
    force_vector = a * particle_body.frame.x
    force_point = particle_body.force_list[0][0]

    frame = particle_body.frame
    assert force_point.vel(frame) == point.vel(frame)
    assert force_point.pos_from(force_point) == point.pos_from(force_point)

    assert particle_body.force_list[0][1] == force_vector


def test_body_add_force():
    # Body with RigidBody.
    rigidbody_masscenter = Point('rigidbody_masscenter')
    rigidbody_mass = Symbol('rigidbody_mass')
    rigidbody_frame = ReferenceFrame('rigidbody_frame')
    body_inertia = inertia(rigidbody_frame, 1, 0, 0)
    rigid_body = Body('rigidbody_body', rigidbody_masscenter, rigidbody_mass,
               rigidbody_frame, body_inertia)

    l = Symbol('l')
    Fa = Symbol('Fa')
    point = rigid_body.masscenter.locatenew(
        'rigidbody_body_point0',
        l * rigid_body.frame.x)
    point.set_vel(rigid_body.frame, 0)
    rigid_body.add_force((0, 0, Fa), point)
    assert len(rigid_body.force_list) == 1
    force_vector = Fa * rigid_body.frame.z
    force_point = rigid_body.force_list[0][0]
    frame = rigid_body.frame
    assert force_point.vel(frame) == point.vel(frame)
    assert force_point.pos_from(force_point) == point.pos_from(force_point)

    assert rigid_body.force_list[0][1] == force_vector
