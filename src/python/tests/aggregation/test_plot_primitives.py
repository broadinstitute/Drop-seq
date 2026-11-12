# MIT License
#
# Copyright 2026 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
The colours of the QC pages are part of their look and must be the same as in the legacy R reports, so these tests read
the fill that each bar or point is drawn with.  The expected colours of the default hue palette were measured in the
PDFs of the R reports (and, for 3 colours, are ggplot2's well-known ones).
"""
import unittest

import pandas as pd

from dropseq.aggregation.plots import plot_primitives as pp

# ggplot2's default hue palette, by number of levels.  7 and 4 are the cell types and libraries of the R baseline.
R_HUE = {
    3: ["#F8766D", "#00BA38", "#619CFF"],
    4: ["#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"],
    7: ["#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7"],
}


def fills(plot, layer):
    """
    :return: dict from the y value of each element of the layer to the colour it is drawn with
    """
    plot.draw()
    data = plot.layers[layer].data
    return {round(y, 6): fill.upper() for y, fill in zip(data["y"], data["fill"])}


class TestColours(unittest.TestCase):
    def test_default_hue_palette_is_assigned_to_the_sorted_labels_whatever_the_bar_order(self):
        for n, expected in R_HUE.items():
            # Alphabetical order is the reverse of the order of the bars, which are sorted by descending value.
            labels = [f"type{i}" for i in range(n)]
            df = pd.DataFrame({"cell_type": labels, "value": [float(i + 1) for i in range(n)]})
            plot = pp.sorted_bar_with_median(df, "cell_type", "value", "# nuclei", None, "ID", "subtitle",
                                             xlabel="Cell type", show_median=False)
            drawn = fills(plot, 0)
            for i, label in enumerate(labels):
                self.assertEqual(drawn[float(i + 1)], expected[i], f"{n} cell types, {label}")

    def test_donors_have_the_donor_palette_in_sorted_order(self):
        donors = ["MS2", "MD1", "UM3", "MD0"]
        df = pd.DataFrame({"donor": donors, "value": [4.0, 3.0, 2.0, 1.0]})
        plot = pp.sorted_bar_with_median(df, "donor", "value", "# nuclei", pp.donor_palette(donors), "ID", "subtitle")
        drawn = fills(plot, 0)
        for value, donor in zip([4.0, 3.0, 2.0, 1.0], donors):
            self.assertEqual(drawn[value], pp.PALETTE[sorted(donors).index(donor)].upper(), donor)

    def test_library_points_have_the_default_hue_palette_whatever_the_row_order(self):
        libraries = [f"lib_rxn{i}" for i in range(1, 5)]
        df = pd.DataFrame({"cell_type": "NN-IMN", "library": libraries[::-1], "value": [4.0, 3.0, 2.0, 1.0]})
        plot = pp.box_with_points(df, "cell_type", "value", "# nuclei", "ID", "subtitle", "Cell type",
                                  color_by="library", legend_title="Library")
        drawn = fills(plot, 1)
        for library, value in zip(libraries[::-1], [4.0, 3.0, 2.0, 1.0]):
            self.assertEqual(drawn[value], R_HUE[4][libraries.index(library)], library)


if __name__ == '__main__':
    unittest.main()
