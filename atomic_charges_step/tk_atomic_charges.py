# -*- coding: utf-8 -*-

"""The graphical part of an Atomic Charges step"""

import logging
import pprint  # noqa: F401
import tkinter as tk

import atomic_charges_step  # noqa: F401
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_widgets as sw  # noqa: F401

logger = logging.getLogger(__name__)


class TkAtomicCharges(seamm.TkNode):
    """The graphical part of an Atomic Charges step in a flowchart.

    See Also
    --------
    AtomicCharges, AtomicChargesParameters
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
    ):
        """Initialize the graphical Tk Atomic Charges node.

        Parameters
        ----------
        tk_flowchart: Tk_Flowchart
            The graphical flowchart that we are in.
        node: Node
            The non-graphical node for this step.
        canvas: tk.Canvas
            The Tk canvas to draw on.
        x, y: float
            The position of the node's center on the canvas.
        w, h: float
            The node's graphical width and height, in pixels.
        """
        self.dialog = None

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
        )

    def create_dialog(self):
        """Create the dialog, building the widgets from the parameters.

        See Also
        --------
        TkAtomicCharges.reset_dialog
        """
        frame = super().create_dialog(title="Atomic Charges")

        # Shortcut for parameters
        P = self.node.parameters

        # Create the widgets for the (non-internal) parameters
        for key in P:
            if key[0] != "_" and key not in ("results",):
                self[key] = P[key].widget(frame)

        # Make the dialog reactive to the controlling choices
        for item in ("method", "density source"):
            w = self[item]
            w.combobox.bind("<<ComboboxSelected>>", self.reset_dialog)
            w.combobox.bind("<Return>", self.reset_dialog)
            w.combobox.bind("<FocusOut>", self.reset_dialog)

        self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Lay out the widgets, showing only those relevant to the choices.

        Parameters
        ----------
        widget : Tk Widget = None
            The widget that triggered the reset, if any.
        """
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        method = self["method"].get()
        source = self["density source"].get()

        row = 0
        self["method"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
        row += 1
        self["charge label"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
        row += 1
        self["density source"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
        row += 1

        # Only show the explicit-files field when not taking it from the
        # previous step.
        if "files" in source:
            self["density files"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            row += 1

        # The reference atomic densities are only used by DDEC6.
        if "DDEC6" in method:
            self["atomic densities directory"].grid(
                row=row, column=0, columnspan=2, sticky=tk.EW
            )
            row += 1

        sw.align_labels(
            [self[k] for k in self.node.parameters if k[0] != "_" and k != "results"],
            sticky=tk.E,
        )

        return row
