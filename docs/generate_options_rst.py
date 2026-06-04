#!/usr/bin/env python
"""
Auto-generate the Options reference page (current_mdv.rst) from VolFe source code.

This script introspects the VolFe module to extract all model parameters,
their defaults, available options, and descriptions, then generates a
comprehensive RST reference page.

Usage:
    python generate_options_rst.py

Kayla Iacovino with function mapping by Claude.
"""

import inspect
import re
import sys
import os

# Ensure the VolFe package is importable
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import VolFe.model_dependent_variables as mdv

# ──────────────────────────────────────────────────────────────────────────────
# Static intro text for the RST page
# ──────────────────────────────────────────────────────────────────────────────

INTRO_RST = r"""=========================
Options
=========================

Various parameters in VolFe can be calculated using different models - these include parameters like solubility functions, equilibrium constants, fugacity coefficients, isotope fractionation factors, etc.
For solubility functions, oxygen fugacity to Fe\ :sup:`3+`/Fe\ :sub:`T` relationships, and sulfur saturation conditions, this is how VolFe includes the effect of melt compositions on calculations.

There are lots of models already available in VolFe, which can be found in the API reference section.
Others can be added as they become available (see :doc:`Add your own <add_your_own>` in the Worked Examples section) - let us know if you have a new model to be added!

Additionally, there are various options for how the calculations are done in VolFe.
This varies from which species are treated as insoluble in the melt to whether a csv of the results is created at the end of the calculation.

How to set options
------------------

Options are passed to calculation functions as a pandas DataFrame created from a list of ``[parameter_name, option_value]`` pairs using :func:`~VolFe.model_dependent_variables.make_df_and_add_model_defaults`:

.. code-block:: python

   import VolFe as vf

   my_models = [
       ['carbon dioxide', 'Basalt_Dixon97'],
       ['gassing_style', 'open'],
   ]
   models = vf.make_df_and_add_model_defaults(my_models)

Any parameters not specified will use the default option. The ``models`` DataFrame is then passed to calculation functions (e.g., ``vf.calc_Pvsat(setup, models=models)``).

Naming conventions
------------------

Here are some hints to the naming convention used for functions, where specific model options for these model dependent variables can be found:

- Equilibrium constants: functions starting with ``K`` (e.g., ``KCOHg()`` is the function for calculating the equilibrium constants for CH\ :sub:`4` + 2O\ :sub:`2` = CO\ :sub:`2` + 2H\ :sub:`2`\ O).

- Fugacity coefficients: functions starting with ``y_`` (e.g., ``y_CH4()`` is the function for calculating the fugacity coefficient for CH\ :sub:`4`).

- Solubility functions: functions starting with ``C_`` (e.g., ``C_CH4()`` is the function for calculating the solubility function for CH\ :sub:`4`).

- Oxygen fugacity and Fe\ :sup:`3+`/Fe\ :sub:`T`: functions are ``FMQ()``, ``NNO()``, ``fO22Fe3FeT()``, and ``f_O2()``.

- Sulfide/sulfate content at sulfide/anhydrite saturation: functions are ``SCSS()`` and ``SCAS()``.

- Isotope fractionation factors: functions are ``alpha_`` (e.g., ``alpha_S_SO2v_S6pm`` is the function for calculating the isotopic fractionation factor between SO\ :sub:`2` in the vapor and S\ :sup:`6+` in the melt).


All available options
---------------------

Below is a comprehensive list of all model parameters, organized by category.
For each parameter, the default option is shown along with all available alternatives.

"""

# ──────────────────────────────────────────────────────────────────────────────
# Category display names and RST heading characters
# ──────────────────────────────────────────────────────────────────────────────

# Map from the raw category strings in the docstring to display names
CATEGORY_DISPLAY = {
    "Specifying species": "Species",
    "Oxygen fugacity": "Oxygen fugacity",
    "Models for solubility and speciation constants": "Solubility constants",
    "Saturation conditions": "Saturation conditions",
    "Fugacity coefficients": "Fugacity coefficients",
    "Equilibrium constants": "Equilibrium constants",
    "Degassing calculation": "Degassing calculation",
    "Other": "Other",
    "In development": "In development",
}

# Map from docstring function references to actual attribute names in mdv module.
# Most map directly but a few differ.
FUNC_NAME_MAP = {
    "fO2": "f_O2",
    "K_COm": "KCOm",
    "K_HOm": "KHOm",
    "KOCg": "KCOg",
}

# Map from docstring parameter names to default_models index names, for cases
# where they differ.
PARAM_NAME_MAP = {
    "KOCg": "KCOg",
}


# ──────────────────────────────────────────────────────────────────────────────
# Parsing helpers
# ──────────────────────────────────────────────────────────────────────────────


def get_defaults():
    """Return dict of parameter_name -> default_value from default_models."""
    return {idx: row["option"] for idx, row in mdv.default_models.iterrows()}


def parse_main_docstring():
    """
    Parse the make_df_and_add_model_defaults docstring.

    Returns a list of (category_name, params) tuples where params is a list of
    dicts with keys: name, description, func_ref, inline_options.
    """
    docstring = inspect.getdoc(mdv.make_df_and_add_model_defaults)
    if docstring is None:
        raise RuntimeError(
            "Could not get docstring from make_df_and_add_model_defaults"
        )

    # Find the "Model Parameters and Options" section
    marker = "Model Parameters and Options"
    idx = docstring.find(marker)
    if idx == -1:
        raise RuntimeError("Could not find 'Model Parameters and Options' in docstring")

    # Skip past the section header and its underline
    section_text = docstring[idx:]
    lines = section_text.split("\n")
    # Skip the header line and the dashes line
    start = 0
    for i, line in enumerate(lines):
        if line.strip().startswith("---"):
            start = i + 1
            break
    lines = lines[start:]

    categories = []
    current_category = None
    current_params = []
    current_param = None

    for line in lines:
        # Check for category header: ### Category Name ###
        cat_match = re.match(r"\s*###\s+(.+?)\s+###", line)
        if cat_match:
            # Save previous category
            if current_category is not None:
                if current_param is not None:
                    current_params.append(current_param)
                    current_param = None
                categories.append((current_category, current_params))
            current_category = cat_match.group(1).strip()
            current_params = []
            continue

        # Check for parameter definition: "param_name: description" or just "param_name"
        # Must start at the beginning of the line (no leading whitespace)
        param_match = re.match(r"^([A-Za-z][\w\s]*?):\s*(.*)", line)
        if param_match and not line.startswith(" ") and not line.startswith("\t"):
            # Save previous param
            if current_param is not None:
                current_params.append(current_param)

            name = param_match.group(1).strip()
            desc = param_match.group(2).strip()
            current_param = {
                "name": name,
                "description": desc,
                "func_ref": None,
                "inline_options": [],
            }
            continue

        # In development parameters
        bare_match = re.match(r"^([a-z_][\w\s]*?)$", line)
        if bare_match and not line.startswith(" ") and current_category:
            if current_param is not None:
                current_params.append(current_param)
            name = bare_match.group(1).strip()
            current_param = {
                "name": name,
                "description": "",
                "func_ref": None,
                "inline_options": [],
            }
            continue

        if current_param is None:
            continue

        # Get reference to function as in doc strings
        func_match = re.search(r"See function (\w+) for options", line)
        if func_match:
            current_param["func_ref"] = func_match.group(1)
            continue

        # Check for inline option: "- 'option_name' [default] description"
        opt_match = re.match(
            r"""\s*-\s+['"]([^'"]+)['"]\s*(?:\[default\])?\s*(.*)""", line
        )
        if opt_match:
            opt_name = opt_match.group(1)
            opt_desc = opt_match.group(2).strip()
            is_default = "[default]" in line
            current_param["inline_options"].append(
                {
                    "name": opt_name,
                    "description": opt_desc,
                    "is_default": is_default,
                }
            )
            continue

        # Check for "default: 'value'" pattern (used in "In development")
        default_match = re.match(r"\s+default:\s+['\"]?(.+?)['\"]?\s*$", line)
        if default_match:
            val = default_match.group(1).strip()
            current_param["inline_options"].append(
                {
                    "name": val,
                    "description": "",
                    "is_default": True,
                }
            )
            continue

        # Append continuation text to description
        stripped = line.strip()
        if (
            stripped
            and current_param["description"]
            and not stripped.startswith("Only one")
        ):
            if not stripped.startswith("For now"):
                current_param["description"] += " " + stripped

    # Save last parameter and category
    if current_param is not None:
        current_params.append(current_param)
    if current_category is not None:
        categories.append((current_category, current_params))

    return categories


def parse_function_options(func_name, param_name=None):
    """
    Parse the "Model options for ..." section from a function's docstring.

    If param_name is given, only return options from the section whose header
    matches that parameter name (some functions like C_H2O have multiple
    "Model options for" sections for different parameters).

    Returns a list of dicts with keys: name, description, is_default.
    """
    # Resolve function name mapping
    actual_name = FUNC_NAME_MAP.get(func_name, func_name)

    func = getattr(mdv, actual_name, None)
    if func is None:
        print(
            f"  Warning: function '{actual_name}' not found in module", file=sys.stderr
        )
        return []

    docstring = inspect.getdoc(func)
    if docstring is None:
        print(f"  Warning: no docstring for '{actual_name}'", file=sys.stderr)
        return []

    # Collect options per section header
    all_sections = {}
    current_section_name = None
    current_options = []
    in_options_section = False

    for line in docstring.split("\n"):
        # Detect start of a "Model options for" section
        header_match = re.match(r"Model options for\s+['\"]?(.+?)['\"]?\s*$", line)
        if header_match:
            # Save previous section if any
            if current_section_name is not None:
                all_sections[current_section_name] = current_options
            current_section_name = header_match.group(1).strip()
            current_options = []
            in_options_section = True
            continue

        # Skip separator lines
        if in_options_section and re.match(r"^-+$", line.strip()):
            continue

        # Parse option lines
        if in_options_section:
            opt_match = re.match(
                r"""\s*-\s+['"]([^'"]+)['"]\s*(?:\[default\])?\s*(.*)""", line
            )
            if opt_match:
                current_options.append(
                    {
                        "name": opt_match.group(1),
                        "description": opt_match.group(2).strip(),
                        "is_default": "[default]" in line,
                    }
                )
            elif line.strip() == "" or line.strip().startswith("Only one"):
                # Blank line or note — could end this section or just be a gap
                pass
            elif not line.startswith(" ") and not line.startswith("-") and line.strip():
                # Non-indented non-option line means we've left the section
                in_options_section = False

    # Save last section
    if current_section_name is not None:
        all_sections[current_section_name] = current_options

    # If a specific param_name was requested, return only that section
    if param_name and param_name in all_sections:
        return all_sections[param_name]

    # Otherwise return all options combined (for functions with a single section)
    all_opts = []
    for opts in all_sections.values():
        all_opts.extend(opts)
    return all_opts


def format_param_rst(param, default_value, func_options=None, is_dev=False):
    """
    Format a single parameter as RST.

    Returns a string of RST content for this parameter.
    """
    name = param["name"]
    desc = param["description"]

    lines = []

    # Parameter header: **name** *(default: value)*
    default_str = str(default_value) if default_value is not None else "?"
    lines.append(f"**{name}** *(default: {default_str})*")

    # For "In development" params, just show the header with default — no
    # description or option bullets needed.
    if is_dev:
        lines.append("")
        return "\n".join(lines)

    # Description
    if desc:
        # Clean up description: remove "See function X for options." since we
        # add a cross-reference separately
        desc_clean = re.sub(r"\s*See function \w+ for options\.?\s*", "", desc).strip()
        # Escape bare asterisks that RST would interpret as emphasis
        desc_clean = re.sub(r"(?<!\S)\*(?!\*)", r"\*", desc_clean)
        if desc_clean:
            func_ref = param.get("func_ref")
            if func_ref:
                actual_func = FUNC_NAME_MAP.get(func_ref, func_ref)
                desc_clean += (
                    f" See :func:`~VolFe.model_dependent_variables.{actual_func}`."
                )
            lines.append(f"   {desc_clean}")
        elif param.get("func_ref"):
            actual_func = FUNC_NAME_MAP.get(param["func_ref"], param["func_ref"])
            lines.append(
                f"   See :func:`~VolFe.model_dependent_variables.{actual_func}`."
            )
    elif param.get("func_ref"):
        actual_func = FUNC_NAME_MAP.get(param["func_ref"], param["func_ref"])
        lines.append(f"   See :func:`~VolFe.model_dependent_variables.{actual_func}`.")

    # Options: prefer function docstring options, fall back to inline
    options = func_options if func_options else param.get("inline_options", [])

    if options:
        lines.append("")
        for opt in options:
            opt_desc = opt["description"]
            if opt_desc:
                lines.append(f"   - ``'{opt['name']}'`` — {opt_desc}")
            else:
                lines.append(f"   - ``'{opt['name']}'``")

    lines.append("")
    return "\n".join(lines)


# ──────────────────────────────────────────────────────────────────────────────
# Main generation
# ──────────────────────────────────────────────────────────────────────────────


def generate_rst():
    """Generate the full RST content for current_mdv.rst."""
    defaults = get_defaults()
    categories = parse_main_docstring()

    parts = [INTRO_RST.lstrip()]

    for cat_raw, params in categories:
        cat_display = CATEGORY_DISPLAY.get(cat_raw, cat_raw)

        # Section heading (using ^ for subsection under "All available options")
        parts.append(cat_display)
        parts.append("^" * len(cat_display))
        parts.append("")

        # Special intro for "In development"
        if "development" in cat_raw.lower():
            parts.append(
                "The following parameters are under active development. "
                "For now, leave them at their default values."
            )
            parts.append("")

        is_dev = "development" in cat_raw.lower()

        for param in params:
            # Resolve parameter name for default lookup
            lookup_name = PARAM_NAME_MAP.get(param["name"], param["name"])
            default_val = defaults.get(lookup_name)

            # Get options from function docstring if referenced
            func_options = None
            if param.get("func_ref"):
                func_options = parse_function_options(
                    param["func_ref"], param_name=param["name"]
                )

            rst = format_param_rst(param, default_val, func_options, is_dev=is_dev)
            parts.append(rst)

    return "\n".join(parts)


def main():
    output_path = os.path.join(os.path.dirname(__file__), "current_mdv.rst")
    rst_content = generate_rst()

    with open(output_path, "w") as f:
        f.write(rst_content)

    print(f"Generated {output_path}")
    print(f"  Total characters: {len(rst_content)}")

    # Count parameters
    param_count = rst_content.count("*(default:")
    print(f"  Parameters documented: {param_count}")


if __name__ == "__main__":
    main()
