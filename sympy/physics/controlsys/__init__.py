"""Controlsys is a package for control systems theory. It currently only
supports linear, time invariant systems via the module 'lti'. You can
create linear control systems in state space or transfer function model
representation, transform these representation into one another,
interconnect systems, evaluate the systems etc.
"""

from .lti import StateSpaceModel, TransferFunctionModel
