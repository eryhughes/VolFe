# Tests that default options are being set correctly

import VolFe as vf


def test_model_default_options():
    options = vf.default_models

    assert options.shape == (78, 1), "Wrong number of default options"

    # check a random sample
    assert options.loc["fO2"]["option"] == "Kress91A"
    assert options.loc["SCSS"]["option"] == "ONeill21hyd"
    assert options.loc["density"]["option"] == "DensityX"
    assert options.loc["output csv"]["option"] == "True"


def test_model_options_update():
    """Tests that model choices update appropriately"""

    options = vf.default_models

    assert options.loc["carbon dioxide"]["option"] == "MORB_Dixon95"
    assert options.loc["hydrogen sulfide"]["option"] == "Basalt_Hughes24"
    assert options.loc["y_S2"]["option"] == "Shi92"

    my_models = [
        ["carbon dioxide", "Basalt_Dixon97"],
        ["hydrogen sulfide", "BasalticAndesite_Hughes24"],
        ["y_S2", "ideal"],
    ]

    new_options = vf.make_df_and_add_model_defaults(my_models)

    assert new_options.loc["carbon dioxide"]["option"] == "Basalt_Dixon97", (
        "carbon dioxide model mismatch"
    )
    assert (
        new_options.loc["hydrogen sulfide"]["option"] == "BasalticAndesite_Hughes24"
    ), "hydrogen sulfide model mismatch"
    assert new_options.loc["y_S2"]["option"] == "ideal", (
        "S2 fugacity coefficient mismatch"
    )
