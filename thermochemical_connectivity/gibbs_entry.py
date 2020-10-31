from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.entries.computed_entries import ComputedStructureEntry

from itertools import combinations
import json
import os
from typing import List

import numpy as np
from scipy.interpolate import interp1d
from monty.json import MontyDecoder, MontyEncoder, MSONable

from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure


__author__ = "Matthew McDermott"
__email__ = "mcdermott@lbl.gov"
__date__ = "October 30, 2020"


with open(os.path.join(os.path.dirname(__file__), "g_els.json")) as f:
    G_ELEMS = json.load(f)
with open(os.path.join(os.path.dirname(__file__), "nist_gas_gf.json")) as f:
    G_GASES = json.load(f)
with open(os.path.join(os.path.dirname(__file__), "compounds.json")) as f:
    G_COMPOUNDS = json.load(f)


class GibbsComputedStructureEntry(ComputedStructureEntry):
    """
    An extension to ComputedStructureEntry which includes the estimated Gibbs
    free energy of formation via a machine-learned model.
    """

    def __init__(
        self,
        structure: Structure,
        formation_enthalpy: float,
        temp: float = 300,
        gibbs_model: str = "SISSO",
        correction: float = 0.0,
        energy_adjustments: list = None,
        parameters: dict = None,
        data: dict = None,
        entry_id: object = None,
    ):
        """
        Args:
            structure (Structure): The pymatgen Structure object of an entry.
            formation_enthalpy (float): Formation enthalpy of the entry, calculated
                using phase diagram construction (eV)
            temp (float): Temperature in Kelvin. If temperature is not selected from
                one of [300, 400, 500, ... 2000 K], then free energies will
                be interpolated. Defaults to 300 K.
            gibbs_model (str): Model for Gibbs Free energy. Currently the default (and
                only supported) option is "SISSO", the descriptor created by Bartel et
                al. (2018).
            correction (float): A correction to be applied to the energy. Defaults to 0
            parameters (dict): An optional dict of parameters associated with
                the entry. Defaults to None.
            data (dict): An optional dict of any additional data associated
                with the entry. Defaults to None.
            entry_id: An optional id to uniquely identify the entry.
        """
        self.structure = structure
        self.formation_enthalpy = formation_enthalpy
        self.temp = temp
        self.interpolated = False

        if self.temp < 300 or self.temp > 2000:
            raise ValueError("Temperature must be selected from range: [300, 2000] K.")

        if self.temp % 100:
            self.interpolated = True

        if gibbs_model.lower() == "sisso":
            gibbs_energy = self.gf_sisso()
        else:
            raise ValueError(
                f"{gibbs_model} not a valid model. Please select from [" f"'SISSO']"
            )

        self.gibbs_model = gibbs_model

        super().__init__(
            structure,
            energy=gibbs_energy,
            correction=correction,
            energy_adjustments=energy_adjustments,
            parameters=parameters,
            data=data,
            entry_id=entry_id,
        )

    def gf_sisso(self) -> float:
        """
        Gibbs Free Energy of formation as calculated by SISSO descriptor from Bartel
        et al. (2018). Units: eV (not normalized)

        WARNING: This descriptor only applies to solids. The implementation here
        attempts to detect and use downloaded NIST-JANAF data for common gases (e.g.
        CO2) where possible.

        Reference: Bartel, C. J., Millican, S. L., Deml, A. M., Rumptz, J. R.,
        Tumas, W., Weimer, A. W., … Holder, A. M. (2018). Physical descriptor for
        the Gibbs energy of inorganic crystalline solids and
        temperature-dependent materials chemistry. Nature Communications, 9(1),
        4168. https://doi.org/10.1038/s41467-018-06682-4

        Returns:
            float: Gibbs free energy of formation (eV)
        """
        comp = self.structure.composition

        if comp.is_element:
            return self.formation_enthalpy

        exp_data = False
        if comp.reduced_formula in G_GASES.keys():
            exp_data = True
            data = G_GASES[comp.reduced_formula]
            factor = comp.get_reduced_formula_and_factor()[1]
        elif comp.reduced_formula in G_COMPOUNDS.keys():
            exp_data = True
            data = G_COMPOUNDS[comp.reduced_formula]
            factor = comp.get_reduced_formula_and_factor()[1]

        if exp_data:
            if self.interpolated:
                g_interp = interp1d([int(t) for t in data.keys()], list(data.values()))
                return g_interp(self.temp) * factor
            else:
                return data[str(self.temp)] * factor

        num_atoms = self.structure.num_sites
        vol_per_atom = self.structure.volume / num_atoms
        reduced_mass = self._reduced_mass()

        return (
            self.formation_enthalpy
            + num_atoms * self._g_delta_sisso(vol_per_atom, reduced_mass, self.temp)
            - self._sum_g_i()
        )

    def _sum_g_i(self) -> float:
        """
        Sum of the stoichiometrically weighted chemical potentials of the elements
        at specified temperature, as acquired from "g_els.json".

        Returns:
             float: sum of weighted chemical potentials [eV]
        """
        elems = self.structure.composition.get_el_amt_dict()

        if self.interpolated:
            sum_g_i = 0
            for elem, amt in elems.items():
                g_interp = interp1d(
                    [float(t) for t in G_ELEMS.keys()],
                    [g_dict[elem] for g_dict in G_ELEMS.values()],
                )
                sum_g_i += amt * g_interp(self.temp)
        else:
            sum_g_i = sum(
                [amt * G_ELEMS[str(self.temp)][elem] for elem, amt in elems.items()]
            )

        return sum_g_i

    def _reduced_mass(self) -> float:
        """
        Reduced mass as calculated via Eq. 6 in Bartel et al. (2018)

        Returns:
            float: reduced mass (amu)
        """
        reduced_comp = self.structure.composition.reduced_composition
        num_elems = len(reduced_comp.elements)
        elem_dict = reduced_comp.get_el_amt_dict()

        denominator = (num_elems - 1) * reduced_comp.num_atoms

        all_pairs = combinations(elem_dict.items(), 2)
        mass_sum = 0

        for pair in all_pairs:
            m_i = Composition(pair[0][0]).weight
            m_j = Composition(pair[1][0]).weight
            alpha_i = pair[0][1]
            alpha_j = pair[1][1]

            mass_sum += (alpha_i + alpha_j) * (m_i * m_j) / (m_i + m_j)

        reduced_mass = (1 / denominator) * mass_sum

        return reduced_mass

    @staticmethod
    def _g_delta_sisso(vol_per_atom, reduced_mass, temp) -> float:
        """
        G^delta as predicted by SISSO-learned descriptor from Eq. (4) in
        Bartel et al. (2018).

        Args:
            vol_per_atom (float): volume per atom [Å^3/atom]
            reduced_mass (float) - reduced mass as calculated with pair-wise sum formula
                [amu]
            temp (float) - Temperature [K]

        Returns:
            float: G^delta
        """

        return (
            (-2.48e-4 * np.log(vol_per_atom) - 8.94e-5 * reduced_mass / vol_per_atom)
            * temp
            + 0.181 * np.log(temp)
            - 0.882
        )

    @classmethod
    def from_pd(
        cls, pd, temp=300, gibbs_model="SISSO"
    ) -> List["GibbsComputedStructureEntry"]:
        """
        Constructor method for initializing a list of GibbsComputedStructureEntry
        objects from an existing T = 0 K phase diagram composed of
        ComputedStructureEntry objects, as acquired from a thermochemical database;
        e.g. The Materials Project.

        Args:
            pd (PhaseDiagram): T = 0 K phase diagram as created in pymatgen. Must
                contain ComputedStructureEntry objects.
            temp (int): Temperature [K] for estimating Gibbs free energy of formation.
            gibbs_model (str): Gibbs model to use; currently the only option is "SISSO".

        Returns:
            [GibbsComputedStructureEntry]: list of new entries which replace the orig.
                entries with inclusion of Gibbs free energy of formation at the
                specified temperature.
        """
        gibbs_entries = []
        for entry in pd.all_entries:
            if (
                entry in pd.el_refs.values()
                or not entry.structure.composition.is_element
            ):
                gibbs_entries.append(
                    cls(
                        entry.structure,
                        formation_enthalpy=pd.get_form_energy(entry),
                        temp=temp,
                        correction=0,
                        gibbs_model=gibbs_model,
                        data=entry.data,
                        entry_id=entry.entry_id,
                    )
                )
        return gibbs_entries

    @classmethod
    def from_entries(
        cls, entries, temp=300, gibbs_model="SISSO"
    ) -> List["GibbsComputedStructureEntry"]:
        """
        Constructor method for initializing GibbsComputedStructureEntry objects from
        T = 0 K ComputedStructureEntry objects, as acquired from a thermochemical
        database e.g. The Materials Project.

        Args:
            entries ([ComputedStructureEntry]): List of ComputedStructureEntry objects,
                as downloaded from The Materials Project API.
            temp (int): Temperature [K] for estimating Gibbs free energy of formation.
            gibbs_model (str): Gibbs model to use; currently the only option is "SISSO".

        Returns:
            [GibbsComputedStructureEntry]: list of new entries which replace the orig.
                entries with inclusion of Gibbs free energy of formation at the
                specified temperature.
        """
        from pymatgen.analysis.phase_diagram import PhaseDiagram

        pd = PhaseDiagram(entries)
        return cls.from_pd(pd, temp, gibbs_model)

    def as_dict(self) -> dict:
        """
        :return: MSONAble dict.
        """
        d = super().as_dict()
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["formation_enthalpy"] = self.formation_enthalpy
        d["temp"] = self.temp
        d["gibbs_model"] = self.gibbs_model
        d["interpolated"] = self.interpolated
        return d

    @classmethod
    def from_dict(cls, d) -> "GibbsComputedStructureEntry":
        """
        :param d: Dict representation.
        :return: GibbsComputedStructureEntry
        """
        dec = MontyDecoder()
        return cls(
            dec.process_decoded(d["structure"]),
            d["formation_enthalpy"],
            d["temp"],
            d["gibbs_model"],
            correction=d["correction"],
            energy_adjustments=[
                dec.process_decoded(e) for e in d.get("energy_adjustments", {})
            ],
            parameters={
                k: dec.process_decoded(v) for k, v in d.get("parameters", {}).items()
            },
            data={k: dec.process_decoded(v) for k, v in d.get("data", {}).items()},
            entry_id=d.get("entry_id", None),
        )

    def __repr__(self):
        output = [
            "GibbsComputedStructureEntry {} - {}".format(
                self.entry_id, self.composition.formula
            ),
            "Gibbs Free Energy (Formation) = {:.4f}".format(self.energy),
        ]
        return "\n".join(output)
