import numpy as np
import os
import json
from pymatgen import Composition, Element
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
import numpy as np

import plotly.express as px

from scipy.spatial import HalfspaceIntersection, ConvexHull

import plotly.graph_objects as go
from pymatgen.util.coord import Simplex


with open(os.path.join(os.path.dirname(__file__), "layouts.json")) as f:
    layouts = json.load(f)


class ChempotMap:
    def __init__(self, pd: PhaseDiagram):
        """
        :param pd: Phase diagram object.
        """
        self.pd = pd
        self.n = pd.dim

    def plot(
        self,
        elements: list = None,
        limits: dict = None,
        comps: list = None,
        comps_mode: str = "mesh",
        comps_colors: list = None,
        label_stable: bool = True,
        default_limit: float = -15.0,
    ):
        """
        Create 3D chemical potential diagram for 3 specified elements and project
        down the stability windows for desired compositions.

        :param elements: Element names.
        :param limits: Dictionary of {Element: [lower_bound, upper_bound]} chemical
        potential limits in calculating chemical potential diagram.
        :param comps: Formulas for compositions to project down and highlight. Can
        also include formulas of phases in the diagram.
        :param comps_mode: Choose from ["mesh","lines","mesh+lines"] to determine
        how the projected compositions are shaded.
        :param comps_colors: List of colors used in shading projected comps.
        :param label_stable: Whether or not to add phase annotations.
        :param default_limit: Default lower limit for chemical potentials of excluded elements.

        :return: Plotly figure.
        """

        lims = np.array([[default_limit, 0]] * self.n)

        if not elements:
            elements = self.pd.elements[:3]
        else:
            elements = [Element(e) for e in elements]

        if not comps_colors:
            comps_colors = px.colors.qualitative.Dark2

        for idx, elem in enumerate(self.pd.elements):
            if limits and elem in limits:
                lims[idx, :] = limits[elem]

        elem_indices = [self.pd.elements.index(e) for e in elements]

        data = self.pd.qhull_data
        hyperplanes = np.insert(
            data, [0], (1 - np.sum(data[:, :-1], axis=1)).reshape(-1, 1), axis=1
        )
        hyperplanes[:, -1] = (
            hyperplanes[:, -1] * -1
        )  # flip to all positive energies (due to defn. of hyperplanes)
        entries = self.pd.qhull_entries

        border_hyperplanes = np.array(([[0] * (self.n + 1)] * (2 * self.n)))

        for idx, limit in enumerate(lims):
            border_hyperplanes[2 * idx, idx] = -1
            border_hyperplanes[2 * idx, -1] = limit[0]
            border_hyperplanes[(2 * idx) + 1, idx] = 1
            border_hyperplanes[(2 * idx) + 1, -1] = limit[1]

        hs_hyperplanes = np.vstack([hyperplanes, border_hyperplanes])

        interior_point = np.average(lims, axis=1).tolist()
        hs_int = HalfspaceIntersection(hs_hyperplanes, np.array(interior_point))

        # organize the boundary points by entry
        domains = {entry: [] for entry in entries}
        for intersection, facet in zip(hs_int.intersections, hs_int.dual_facets):
            for v in facet:
                if v < len(entries):
                    this_entry = entries[v]
                    domains[this_entry].append(intersection)

        domains = {k: v for k, v in domains.items() if v}
        domain_vertices = {}
        annotations = []
        font_dict = {"color": "black", "size": 16.0}
        opacity = 0.7

        extra_domains = {}

        for entry, points in domains.items():
            points = np.array(points)
            points_3d = np.array(points[:, elem_indices])
            contains_target_elems = set(entry.composition.elements).issubset(elements)

            if comps:
                if entry.composition.reduced_composition in [
                    Composition(comp).reduced_composition for comp in comps
                ]:
                    domains[entry] = None
                    extra_domains[entry] = points_3d

                    if contains_target_elems:
                        domains[entry] = points_3d
                    else:
                        continue

            if not contains_target_elems:
                domains[entry] = None
                continue

            try:
                domain = ConvexHull(points_3d)
                ann_loc = np.mean(points_3d.T, axis=1)
            except:
                points_2d, v, w = self.simple_pca(points_3d, k=2)
                domain = ConvexHull(points_2d)
                centroid_2d = self.get_centroid_2d(points_2d[domain.vertices])
                ann_loc = centroid_2d @ w.T + np.mean(
                    points_3d.T, axis=1
                )  # recover orig 3D coords from eigenvectors

            simplices = [Simplex(points_3d[indices]) for indices in domain.simplices]

            formula = entry.composition.reduced_formula
            if hasattr(entry, "original_entry"):
                formula = entry.original_entry.composition.reduced_formula

            clean_formula = PDPlotter._htmlize_formula(formula)
            annotation = layouts["default_chempot_annotation_layout"].copy()

            annotation.update(
                {
                    "x": ann_loc[0],
                    "y": ann_loc[1],
                    "z": ann_loc[2],
                    "font": font_dict,
                    "text": clean_formula,
                    "opacity": opacity,
                }
            )
            annotations.append(annotation)
            domains[entry] = simplices
            domain_vertices[entry] = points_3d

        x, y, z = [], [], []

        for phase, simplexes in domains.items():
            if simplexes:
                for s in simplexes:
                    x.extend(s.coords[:, 0].tolist() + [None])
                    y.extend(s.coords[:, 1].tolist() + [None])
                    z.extend(s.coords[:, 2].tolist() + [None])

        layout = layouts["default_chempot_layout_3d"].copy()
        layout["scene"].update(
            {
                "xaxis": self.get_chempot_axis_layout(elements[0]),
                "yaxis": self.get_chempot_axis_layout(elements[1]),
                "zaxis": self.get_chempot_axis_layout(elements[2]),
            }
        )
        if label_stable:
            layout["scene"].update({"annotations": annotations})

        lines = [
            go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="lines",
                line=dict(color="black", width=4.5),
                showlegend=False,
            )
        ]
        extra_phases = []

        for idx, (entry, coords) in enumerate(extra_domains.items()):
            points_3d = coords[:, :3]
            if "mesh" in comps_mode:
                extra_phases.append(
                    go.Mesh3d(
                        x=points_3d[:, 0],
                        y=points_3d[:, 1],
                        z=points_3d[:, 2],
                        alphahull=0,
                        showlegend=True,
                        lighting=dict(fresnel=1.0),
                        color=comps_colors[idx],
                        name=f"{entry.composition.reduced_formula} (mesh)",
                        opacity=0.13,
                    )
                )
            if "lines" in comps_mode:
                points_2d = points_3d[:, 0:2]
                domain = ConvexHull(points_2d)
                simplexes = [
                    Simplex(points_3d[indices]) for indices in domain.simplices
                ]
                x, y, z = [], [], []
                for s in simplexes:
                    x.extend(s.coords[:, 0].tolist() + [None])
                    y.extend(s.coords[:, 1].tolist() + [None])
                    z.extend(s.coords[:, 2].tolist() + [None])

                extra_phases.append(
                    go.Scatter3d(
                        x=x,
                        y=y,
                        z=z,
                        mode="lines",
                        line={"width": 8, "color": comps_colors[idx]},
                        opacity=1.0,
                        name=f"{entry.composition.reduced_formula} (lines)",
                    )
                )

        layout["scene_camera"] = dict(
            eye=dict(x=0, y=0, z=2.0), projection=dict(type="orthographic")
        )
        fig = go.Figure(lines + extra_phases, layout)

        return fig

    @staticmethod
    def get_chempot_axis_layout(element):
        return dict(
            title=f"μ<sub>{str(element)}</sub> - μ<sub>"
            f"{str(element)}</sub><sup>o</sup> (eV)",
            titlefont={"size": 30},
            gridcolor="#e8e8e8",
            gridwidth=3.5,
            tickfont={"size": 16},
            ticks="inside",
            ticklen=14,
            showline=True,
            backgroundcolor="rgba(0,0,0,0)",
        )

    @staticmethod
    def simple_pca(data, k=2):
        """ Simple implementation of Principal Component Analaysis"""
        data = data - np.mean(data.T, axis=1)  # centering the exp_data
        cov = np.cov(data.T)  # calculating covariance matrix
        v, w = np.linalg.eig(cov)  # performing eigendecomposition
        idx = v.argsort()[::-1]  # sorting the components
        v = v[idx]
        w = w[:, idx]
        scores = data.dot(w[:, :k])

        return scores, v[:k], w[:, :k]

    @staticmethod
    def get_centroid_2d(vertices):
        """ Simple implementation of 2D centroid formula for annotations. Vertices must
        be provided in order going around around the polygon"""
        n = len(vertices)
        cx = 0
        cy = 0
        a = 0

        for i in range(0, n - 1):
            xi = vertices[i, 0]
            yi = vertices[i, 1]
            xi_p = vertices[i + 1, 0]
            yi_p = vertices[i + 1, 1]
            common_term = xi * yi_p - xi_p * yi

            cx += (xi + xi_p) * common_term
            cy += (yi + yi_p) * common_term
            a += common_term

        prefactor = 0.5 / (6 * a)
        return np.array([prefactor * cx, prefactor * cy])
