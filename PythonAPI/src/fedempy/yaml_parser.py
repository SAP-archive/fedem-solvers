# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

"""
This module provides functionality for creating a Fedem model based on a
YAML-formatted input file. The following keywords are recognized when parsing
the input file:

| *target_file* - specifies the name of the fmm-file to be created
| *source_file* - specifies an fmm-file to be used as template for the model.
| *file_exists* - specified what to do if the `target_file`already exists:
|                 USE_IT - opens this file for editing existing model
|                 OVERWRITE! - ignores and overwrites any existing file
|                 STOP! - aborts the execution
|                 index - creates a new non-existing file name
| *fe_parts* - imports specified FE data files as a Parts
| *triads* - creates triads at specified global points
| *triads_from_fe_parts* - creates triads at nodal points on specified FE Parts
| *beam_materials* - creates Beam material objects
| *beam_sections* - creates Beam cross section objects
| *beams* - creates Beam element objects
| *joints* - creates Joint objects
| *loads* - create Load objects
| *virtual_sensors* - creates Sensor objects for extraction of response data
| *functions* - creates Function objects (including input functions)
| *spring_dampers* - creates Axial spring/damper objects
| *edit_fe_parts* - modifies the properties of FE Part objects
| *edit_triads* - modifies the properties of Triad objects
| *edit_joints* - modifies the properties of Joint objects
| *edit_settings* - modifies the solver settings

Each keyword (except the first three) are followed by an arbitrary number of
lines, where each line defines one object to be created (or edited).
The syntax of the input can be either list-based or dictionary based.

This module can also be launched directly using the syntax

| ``python -m fedempy.yaml_parser -f mymodel.yaml``

It will then invoke the method :meth:`yaml_parser.ModelYAML.main`
on the specified input file (`mymodel.yaml`).

"""

from argparse import ArgumentParser
from collections import OrderedDict
from json import dump
from os import getcwd, makedirs, path
from pathlib import Path
from shutil import copy2

from yaml import FullLoader, load

from fedempy.enums import FmDof, FmDofStat, FmLoadType, FmType, FmVar
from fedempy.fmm_solver import FedemException, FmmSolver
from fedempy.modeler import FedemModeler


def _check_target(rel_path, fname, exists):
    """
    Check for existing target file.
    """
    if path.exists(rel_path + fname):
        # Manage already existing target model file
        if exists == "STOP!":
            raise FedemException(f"Target file {fname} alread exists!")
        if exists == "OVERWRITE!":
            print(f"  ** Existing file {fname} will be overwritten.")
        elif exists == "USE_IT":
            print(f"   * Updating existing model file {fname}.")
        elif exists == "index":
            for i in range(1000):
                new_name = f"{fname[:-4]}-i{i:03}{fname[-4:]}"
                if not path.exists(rel_path + new_name):
                    break
            print(f"  ** Target file {fname} exists, using {new_name}.")
            fname = new_name
        else:
            raise FedemException(f"Invalid command ({exists})")

    return rel_path + fname


def _merge(target, source):
    """
    Merges the `source` dictionary into the `target` dictionary.
    Entries for already existing keys shall _not_ be overwritten.
    Therefore, we can't use the update() method here.
    """
    for _key, _val in source.items():
        # Account for that "triads" and "triads_from_fe_parts" are equivalent
        if _key in target:
            my_key = _key
        elif _key == "triads" and "triads_from_fe_parts" in target:
            my_key = "triads_from_fe_parts"
        elif _key == "triads_from_fe_parts" and "triads" in target:
            my_key = "triads"
        else:
            target[_key] = _val
            continue

        for __key, __val in _val.items():
            if __key not in target[my_key]:
                target[my_key][__key] = __val

    return target


def _process(generator, properties):
    """
    Generic function for invoking the method `generator`
    on the data dictionary `properties`.
    """
    for name, data in properties.items():
        if name.startswith("-"):
            print("   * Comment line:", name, data)
        else:
            print(f"\t{name}")
            if generator(name, data) <= 0:
                raise FedemException("Model generation failure.")


class ModelYAML:
    """
    This class contains methods for parsing a YAML-formated input file
    and creating a Fedem model from it. It is an extension of the class
    :class:`modeler.FedemModeler`.

    Parameters
    ----------
    input_file : str
        Absolute path to the YAML-formated input file

    Methods
    -------
    build:
        Loads and parses all model attributes into the model instance
    save:
        Saves the model to the input file specified FEDEM model file
    solve:
        Executes the dynamics solver on the created model.
    """

    def __init__(self, input_file):
        """
        Constructor.
        """
        self.model_dir = Path(input_file).parent.resolve()
        self.model_file = None
        self.model_input = {}

        print("   * Parsing YAML input file", input_file)
        with open(input_file) as yml_file:
            self.model_input = load(yml_file, Loader=FullLoader)

        if self.model_input:
            # Check for inclusion of other yaml-file(s) into parent file
            sub_models = self._include(self.model_input.pop("include", {}))
            if sub_models:
                # Merge the content of the included files with the current
                # yaml-file content. Notice that any existing values in
                # the current dictionary will not be overwritten.
                self.model_input = _merge(self.model_input, sub_models)

            # Extract some file management keywords
            rel_path = self.model_input.pop("relative_path", "./")
            def_file = input_file[: input_file.rfind(".")] + ".fmm"
            target_file = self.model_input.pop("target_file", def_file)
            source_file = self.model_input.pop("source_file", None)
            file_exists = self.model_input.pop("file_exists", "index")
            if target_file.startswith("/"):  # Absolute path
                rel_path = ""
            elif getcwd() == "/":
                # Necessary for running in docker (start directory is then "/")
                rel_path = path.join(self.model_dir, rel_path)
            if source_file == rel_path + target_file:
                file_exists = "USE_IT"  # We are amending an existing model
            # Assign the resulting target_file to the [self.model_file] variable
            self.model_file = _check_target(rel_path, target_file, file_exists)
            if rel_path not in ("", "./"):
                makedirs(rel_path, exist_ok=True)

            if (
                source_file is not None
                and path.exists(source_file)
                and file_exists != "USE IT"
            ):
                copy2(source_file, self.model_file)

        # Instanciate a new Fedem model object
        self.model = FedemModeler(self.model_file, file_exists == "OVERWRITE!")

        # Establish keyword-to-method mapping
        self.model_generators = {
            # Make standard FEDEM attributes
            "fe_parts": self._make_fe_part,
            "triads": self._make_triad,
            "triads_from_fe_parts": self._make_triad,
            "beam_materials": self._make_beam_material,
            "beam_sections": self._make_beam_section,
            "beams": self._make_beam,
            "joints": self._make_joint,
            "loads": self._make_load,
            "virtual_sensors": self._make_sensor,
            "functions": self._make_function,
            "spring_dampers": self._make_spring_damper,
            # Edit standard FEDEM attributes
            "edit_fe_parts": self._edit_part,
            "edit_triads": self._edit_triad,
            "edit_joints": self._edit_joint,
            # Edit solver setup
            "edit_settings": self._edit_settings,
        }

    def build(self, dump_file=None):
        """
        Loads and parses all model attributes into the model instance.

        Parameters
        ----------
        dump_file : str, default=None
            Absolute path to dump the model data dictionary to
        """
        if dump_file:
            with open(dump_file, "w") as d_file:
                dump(OrderedDict(sorted(self.model_input.items())), d_file, indent=2)

        model_edit = {}
        for method, properties in self.model_input.items():
            if method.startswith("-"):
                print("   * Comment line:", method)
            elif method not in self.model_generators:
                print("  ** Ignoring unknown keyword:", method)
            elif method.startswith("edit_"):
                # All editing commands will be processed at the end
                model_edit[method] = properties
            else:
                print(f"\n********** Processing {method} **************")
                _process(self.model_generators[method], properties)

        # Now process the edit-commands
        for method, properties in model_edit.items():
            print(f"\n********** Processing {method} **************")
            _process(self.model_generators[method], properties)

    def save(self, save_as=False):
        """
        Saves the model with the defined file name,
        or save it as a new file if `save_as` is set.

        Parameters
        ----------
        save_as : bool, default=False
            If True, the model is stored with a new file name
        """
        new_name = input("\n   * Enter new model name : ") if save_as else None
        if not self.model.save(new_name):
            raise FedemException("Failure saving the Fedem model.")

    def solve(self):
        """
        Executes the dynamics solver on the created model.
        """
        solver = FmmSolver()
        if solver.solve_all(self.model_file, True, True) < 0:
            raise FedemException("Failure solving the Fedem model.")

    def _include(self, yaml_files, level=0):
        """
        Handles the inclusion of other yaml-files.
        """
        if not yaml_files:
            return {}  # No included files

        # Check argument type, only string or a list if strings are accepted
        if isinstance(yaml_files, str):
            yaml_files = [yaml_files]
        elif not isinstance(yaml_files, list):
            raise FedemException("Invalid specification of include files")

        model_input = {}
        for included_file in yaml_files:
            # Parse the included file
            yaml_name = path.join(self.model_dir, included_file)
            print("   *", (" " * level), "Parsing included file", yaml_name)
            submod = load(open(yaml_name), Loader=FullLoader)
            # Recursively check for sub-models
            submod = _merge(submod, self._include(submod.pop("include", {}), level + 1))
            # Merge the sub-models.
            # Notice that any existing values in model_input will not be overwritten.
            model_input = _merge(model_input, submod)

        return model_input

    def _get_triad_id(self, tags, get_all=False):
        """
        Returns the base id of Triad(s) with matching tag.
        """
        triads = self.model.fm_get_objects(FmType.TRIAD, tags)
        if get_all:
            return triads
        if len(triads) > 0:
            return triads[0]

        print(" *** No triad matches the tag(s)", tags)
        return 0

    def _get_joint_id(self, tags, get_all=False):
        """
        Returns the base id of Joint(s) with matching tag.
        """
        joints = self.model.fm_get_objects(FmType.JOINT, tags)
        if get_all:
            return joints
        if len(joints) > 0:
            return joints[0]

        print(" *** No joint matches the tag(s)", tags)
        return 0

    def _get_part_id(self, tags, get_all=False):
        """
        Returns the base id of FE part(s) with matching tag.
        """
        parts = self.model.fm_get_objects(FmType.FEPART, tags)
        if get_all:
            return parts
        if len(parts) > 0:
            return parts[0]

        print(" *** No FE part matches the tag(s)", tags)
        return 0

    def __make_triad(self, name, *vargs, **kwargs):
        """
        Appends a triad to the model.
        """
        data = {}
        for i, key in enumerate(["pos", "rot", "node", "on_part"]):
            if i < len(vargs):
                data[key] = vargs[i]
        data.update(kwargs)

        on_part = self._get_part_id(data["on_part"]) if "on_part" in data else 0

        if "pos" in data:
            pos = data.get("pos", None)
            rot = data.get("rot", None)
            return self.model.make_triad(name, pos, rot, on_part=on_part)

        if "node" in data:
            return self.model.make_triad(name, node=data["node"], on_part=on_part)

        print(" *** Neither position data nor node index is given.")
        return 0

    def _make_triad(self, key, values):
        """
        Appends a triad to the model.
        """
        if isinstance(values, list):
            name = values.pop(0) if isinstance(values[0], str) else ""
            if len(values) == 3 and values[1] == "node" and isinstance(values[2], int):
                part_id = self._get_part_id(values[0])
                base_id = self.model.make_triad(name, node=values[2], on_part=part_id)
            else:
                base_id = self.__make_triad(name, *values)
        elif isinstance(values, dict):
            base_id = self.__make_triad(**values)
        else:
            base_id = 0

        if base_id <= 0:
            print(f" *** Triad {key} with values {values} could not be created.")
        else:
            self.model.fm_tag_object(base_id, key)

        return base_id

    def _make_beam_material(self, key, values):
        """
        Appends a beam material object to the model.
        """
        if isinstance(values, list):
            name = values.pop(0) if isinstance(values[0], str) else ""
            prop = values
        elif isinstance(values, dict):
            name = values.get("name", "")
            prop = [
                values.get("rho", 7850),
                values.get("E", 2.1e11),
                values.get("nu", 0.295),
            ]
        else:
            print(f" *** Material {key} with values {values} has wrong syntax.")
            return 0

        base_id = self.model.make_beam_material(name, prop, key)
        if base_id <= 0:
            print(f" *** Material {key} with values {values} could not be created.")

        return base_id

    def _make_beam_section(self, key, values):
        """
        Appends a beam cross section object to the model.
        """
        if isinstance(values, list):
            name = values.pop(0) if isinstance(values[0], str) else ""
            material_key = values.pop(0)
            parameters = values
        elif isinstance(values, dict):
            name = values.get("name", "")
            material_key = values.get("material", 0)
            parameters = values.get("parameters", {})
        else:
            print(f" *** Beam section {key} with values {values} has wrong syntax.")
            return 0

        material_id = self.model.fm_get_objects(FmType.MAT_PROP, material_key)
        if len(material_id) == 0:
            print(f" *** No beam material {material_key} in the model.")
            return 0

        base_id = self.model.make_beam_section(name, material_id[0], parameters, key)
        if base_id <= 0:
            print(f" *** Beam section {key} with values {values} could not be created.")

        return base_id

    def _make_beam(self, key, values):
        """
        Append a beam element to the model.
        """
        if isinstance(values, list):
            name = values.pop(0) if isinstance(values[0], str) else ""
            triads = self._get_triad_id(values[0], True)
            beam_section = self.model.fm_get_objects(FmType.BEAM_PROP, values[1])
            if len(beam_section) == 0:
                print(f" *** No beam section {values[1]} in the model.")
                return 0
        else:
            print(f" *** Beam {key} with values {values} has wrong syntax")
            return 0

        base_id = self.model.make_beam(name, triads, beam_section[0], key)
        if len(base_id) == 0:
            print(f" *** Beam(s) {key} with values {values} could not be created.")

        return base_id[0]

    def _make_joint(self, key, values):
        """
        Appends a joint to the model.
        """
        name = values.pop(0) if isinstance(values[0], str) else ""
        edit_data = values.pop(-1) if isinstance(values[-1], dict) else False
        values[0] = getattr(FmType, values[0])
        values[1:] = self._get_triad_id(values[1:], True)
        if len(values) < 3:
            values.append(None)
        if len(values) < 4:
            values.append(key)
        else:
            values[3] = key

        base_id = self.model.make_joint(name, *values)
        if base_id <= 0:
            print(f" *** Joint {key} with values {values} could not be created.")
        elif edit_data:
            self._edit_joint(base_id, edit_data)

        return base_id

    def _make_sensor(self, key, values):
        """
        Appends a sensor object to the model.
        """
        name = values.pop(0) if isinstance(values[0], str) else ""
        base_id = self.model.fm_get_objects(FmType.ALL, values[0])
        if len(base_id) < 1:
            print(f" *** Object {values[0]} not found.")
            return 0

        values[0] = base_id[0]
        values[1] = getattr(FmVar, values[1])
        if len(values) > 2:
            values[2] = getattr(FmDof, values[2])
        user_id = self.model.make_sensor(name, *values, tag=key)
        if user_id <= 0:
            print(f" *** Sensor {key} with values {values} could not be created.")

        return user_id

    def _make_load(self, key, values):
        """
        Appends a load to the model.
        """
        name = values.pop(0) if isinstance(values[0], str) else ""
        values[0] = getattr(FmLoadType, values[0])
        values[1] = self._get_triad_id(values[1])
        base_id = self.model.make_load(name, *values, tag=key)
        if base_id <= 0:
            print(f" *** Load {key} with values {values} could not be created.")

        return base_id

    def __make_fe_part(self, key, *vargs, **kwargs):
        """
        Appends an FE part to the model.
        """
        data = {}
        for i, _key in enumerate(
            ["file_name", "name", "Tx", "Ty", "Tz", "Rx", "Ry", "Rz"]
        ):
            if i < len(vargs):
                data[_key] = vargs[i]
        data.update(kwargs)

        library_path = data.pop("library_path", self.model_dir)
        file_name = library_path / data.pop("file_name")
        base_id = self.model.make_fe_part(str(file_name), data.get("name"), tag=key)
        if base_id <= 0:
            print(f" *** FE Part {key} with name {file_name} could not be created.")
            return base_id

        triads = data.pop("triads", {})
        if len(triads) > 0:
            for _key, _kwargs in triads.items():
                self._make_triad(f"{_key}@{key}", **_kwargs, on_part=base_id)

        if len(data) > 0:
            self.model.edit_part(base_id, **data)

        return base_id

    def _make_fe_part(self, key, values):
        """
        Appends an FE part to the model.
        """
        if isinstance(values, list):
            return self.__make_fe_part(key, *values)
        if isinstance(values, dict):
            return self.__make_fe_part(key, **values)

        print(f" *** FE part {key} with values {values} has wrong syntax.")
        return 0

    def _make_spring_damper(self, key, values):
        """
        Appends an axial spring/damper object to the model.
        """
        name = values.pop(0) if isinstance(values[0], str) else ""
        triads = self.model.fm_get_objects(FmType.TRIAD, values[0:2])
        if len(triads) < 2:
            print(f" *** Too few triads ({len(triads)}) for spring/damper {key}")
            return 0

        triads = (triads[0], triads[1])
        bid = self.model.make_spring(name, triads, init_Stiff_Coeff=values[2], tag=key)
        if bid <= 0:
            print(f" *** Spring {key} with values {values} could not be created.")

        if bid > 0 and len(values) > 3:
            self.model.make_damper(name, triads, init_Damp_Coeff=values[3], tag=key)

        return bid

    def _make_function(self, key, values):
        """
        Appends a function object to the model.
        """
        if isinstance(values, str):
            return self.model.make_function(values, tag=key)

        return self.model.make_function(key, **values, tag=key)

    def _edit_triad(self, key, data):
        """
        Modifies the position or properties of a triad.
        """
        if self.model.edit_triad(self._get_triad_id(key, True), **data):
            return 1

        print(f" *** Failed to edit Triad {key} with data {data}.")
        return 0

    def _edit_joint(self, key, data):
        """
        Modifies the position or properties of a joint.
        """
        if "constraints" not in data:
            data["constraints"] = {}
        if "spring" in data:
            for k in data["spring"].keys():
                data["constraints"][k] = FmDofStat.SPRING_DYN
        if "damper" in data:
            for k in data["damper"].keys():
                data["constraints"][k] = FmDofStat.SPRING_DYN

        if isinstance(key, int):
            if self.model.edit_joint(key, **data):
                return 1
        elif self.model.edit_joint(self._get_joint_id(key, True), **data):
            return 1

        print(f" *** Failed to edit Joint {key} with data {data}.")
        return 0

    def _edit_part(self, key, data):
        """
        Modifies the position or properties of an FE part.
        """
        if self.model.edit_part(self._get_part_id(key, True), **data):
            return 1

        print(f" *** Failed to edit Part {key} with data {data}.")
        return 0

    def _edit_settings(self, key, data):
        """
        Setting solver parameters like t_start, t_inc, t_end, t_quasi, e_inc, n_modes.
        Usage in yaml file, e.g.:
        edit_settings:
            solver: { t_start: 0., t_inc: 0.05, t_end: 20., t_quasi: -20., e_inc: 0., n_modes: 0 }
        """
        if key == "solver":
            self.model.solver_setup(**data)
        else:
            print(f"  ** Key {key} with {data} is not implemented yet (ignored).")

        return 1


def main(input_file, dump_file=None, solve=False):
    """
    Main driver.

    Parameters
    ----------
    input_file : str
        Absolute path to the YAML-formated input file
    dump_file : str, default=None
        Absolute path to json-file for dumping the model data dictionary
    solve : bool, default=False
        If True, the dynamics solver is launched on the created model
    """

    model = ModelYAML(input_file)

    model.build(dump_file)
    model.save()
    model.model.close()

    if solve:
        model.solve()


if __name__ == "__main__":
    parser = ArgumentParser(description="PythonAPI YAML parser")
    parser.add_argument("-f", "--input-file", required=True, help="YAML input file")
    parser.add_argument("-d", "--dump-file", help="JSON output file for model dump")
    parser.add_argument("-s", "--solve", action="store_true", help="Solves the model")
    main(**vars(parser.parse_args()))
