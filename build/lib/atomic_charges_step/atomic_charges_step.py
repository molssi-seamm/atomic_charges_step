# -*- coding: utf-8 -*-

import atomic_charges_step


class AtomicChargesStep(object):
    """Helper class needed for the stevedore integration.

    This must provide a `description()` method that returns a dict containing a
    description of this node, and `create_node()` and `create_tk_node()` methods
    for creating the graphical and non-graphical nodes.

    Attributes
    ----------
    my_description : {str, str}
        A human-readable description of this step, with the keys "description",
        "group", and "name".
    """

    my_description = {
        "description": (
            "Compute atomic (partial) charges from the electron density "
            "produced by a preceding quantum-chemistry step."
        ),
        "group": "Analysis",
        "name": "Atomic Charges",
    }

    def __init__(self, flowchart=None, gui=None):
        """Initialize this helper class, which is used by
        the application via stevedore to get information about
        and create node objects for the flowchart
        """
        pass

    def description(self):
        """Return a description of what this step does.

        Returns
        -------
        description : dict(str, str)
        """
        return AtomicChargesStep.my_description

    def create_node(self, flowchart=None, **kwargs):
        """Create and return the new node object.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart this node is part of.

        Returns
        -------
        AtomicCharges
        """
        return atomic_charges_step.AtomicCharges(flowchart=flowchart, **kwargs)

    def create_tk_node(self, canvas=None, **kwargs):
        """Create and return the graphical Tk node object.

        Parameters
        ----------
        canvas : tk.Canvas
            The Tk Canvas widget

        Returns
        -------
        TkAtomicCharges
        """
        return atomic_charges_step.TkAtomicCharges(canvas=canvas, **kwargs)
