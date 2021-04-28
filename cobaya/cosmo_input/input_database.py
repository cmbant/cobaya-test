# Notes about writing more presets:
# ---------------------------------
# - all parameter names below are PLANCK parameter names. They are substituted by the
#   theory-code-specific ones in `create_input`
# - don't use extra_args for precision parameters! because if the same precision param
#   is mentioned twice at the same time in different fields with different values, there
#   is no facility to take the max (or min). Instead, codify precision needs in terms of
#   requirements in the .must_provide method of the cosmo code.

# Global
from copy import deepcopy

# Local
from cobaya.conventions import kinds, _params
from cobaya.typing import InfoDict

_camb = "camb"
_classy = "classy"
_desc = "desc"
_comment = "note"
_extra_args = "extra_args"
_error_msg = "error_msg"
_none = "(None)"

# Theory codes
theory: InfoDict = {_camb: None, _classy: None}

# Primordial perturbations
primordial: InfoDict = dict([
    ("SFSR", {
        _desc: "Adiabatic scalar perturbations, power law spectrum",
        kinds.theory: {_camb: None, _classy: None},
        _params: dict([
            ["logA", dict([
                ("prior", dict([("min", 1.61), ("max", 3.91)])),
                ("ref", dict([("dist", "norm"), ("loc", 3.05), ("scale", 0.001)])),
                ("proposal", 0.001), ("latex", r"\log(10^{10} A_\mathrm{s})"),
                ("drop", True)])],
            ("As", dict([
                ("value", "lambda logA: 1e-10*np.exp(logA)"),
                ("latex", r"A_\mathrm{s}")])),
            ("ns", dict([
                ("prior", dict([("min", 0.8), ("max", 1.2)])),
                ("ref", dict([("dist", "norm"), ("loc", 0.965), ("scale", 0.004)])),
                ("proposal", 0.002), ("latex", r"n_\mathrm{s}")]))])})])
primordial.update(dict([
    ["SFSR_DESpriors", {
        _desc: "Adiabatic scalar perturbations, power law + running spectrum "
               "-- DESpriors",
        kinds.theory: {_camb: None, _classy: None},
        _params: dict([
            ["As_1e9", dict([
                ["prior", dict([["min", 0.5], ["max", 5]])],
                ["ref",
                 dict([["dist", "norm"], ["loc", 2.1], ["scale", 0.5]])],
                ["proposal", 0.25], ["latex", r"10^9 A_\mathrm{s})"],
                ["drop", True], ["renames", "A"]])],
            ["As", dict([
                ["value", "lambda As_1e9: 1e-9 * As_1e9"],
                ["latex", r"A_\mathrm{s}"]])],
            ["ns", primordial["SFSR"][_params]["ns"]]])}]]))
primordial.update(dict([
    ["SFSR_run", {
        _desc: "Adiabatic scalar perturbations, power law + running spectrum",
        kinds.theory: {_camb: None, _classy: None},
        _params: dict(
            (list(primordial["SFSR"][_params].items()) +
             [["nrun", dict([
                 ["prior", dict([["min", -1], ["max", 1]])],
                 ["ref",
                  dict([["dist", "norm"], ["loc", 0], ["scale", 0.005]])],
                 ["proposal", 0.001], ["latex", r"n_\mathrm{run}"]])]]))}]]))
primordial.update(dict([
    ["SFSR_runrun", {
        _desc: "Adiabatic scalar perturbations, power law + 2nd-order running spectrum",
        kinds.theory: {_camb: None, _classy: None},
        _params: dict(
            (list(primordial["SFSR_run"][_params].items()) +
             [["nrunrun", dict([
                 ["prior", dict([["min", -1], ["max", 1]])],
                 ["ref",
                  dict([["dist", "norm"], ["loc", 0], ["scale", 0.002]])],
                 ["proposal", 0.001],
                 ["latex", r"n_\mathrm{run,run}"]])]]))}]]))
primordial.update(dict([
    ["SFSR_t", {
        _desc: "Adiabatic scalar+tensor perturbations, "
               "power law spectrum (inflation consistency)",
        kinds.theory: {_camb: {_extra_args: {"nt": None}},
                       _classy: {_extra_args: {"n_t": "scc", "alpha_t": "scc"}}},
        _params: dict(
            (list(primordial["SFSR"][_params].items()) +
             [["r", dict([
                 ["prior", dict([["min", 0], ["max", 3]])],
                 ["ref",
                  dict([["dist", "norm"], ["loc", 0], ["scale", 0.03]])],
                 ["proposal", 0.03], ["latex", r"r_{0.05}"]])]]))}]]))
primordial.update(dict([
    ["SFSR_t_nrun", {
        _desc: "Adiabatic scalar+tensor perturbations, "
               "power law + running spectrum (inflation consistency)",
        kinds.theory: {_camb:
                           {_extra_args: primordial["SFSR_t"][kinds.theory][_camb][
                               _extra_args]},
                       _classy:
                           {_extra_args: primordial["SFSR_t"][kinds.theory][_classy][
                               _extra_args]}},
        _params: dict(
            (list(primordial["SFSR_run"][_params].items()) +
             list(primordial["SFSR_t"][_params].items())))}]]))

# Geometry
geometry: InfoDict = dict([
    ["flat", {
        _desc: "Flat FLRW universe",
        kinds.theory: {_camb: None, _classy: None}}],
    ["omegak", {
        _desc: "FLRW model with varying curvature (prior on Omega_k)",
        kinds.theory: {_camb: None, _classy: None},
        _params: dict([
            ["omegak", dict([
                ["prior", dict([["min", -0.3], ["max", 0.3]])],
                ["ref",
                 dict([["dist", "norm"], ["loc", -0.009], ["scale", 0.001]])],
                ["proposal", 0.001], ["latex", r"\Omega_k"]])]])}], ])

# Hubble parameter constraints
H0_min, H0_max = 20, 100
hubble: InfoDict = dict([
    ["H", {
        _desc: "Hubble parameter",
        kinds.theory: {_camb: None, _classy: None},
        _params: dict([
            ["H0", dict([
                ["prior", dict([["min", H0_min], ["max", H0_max]])],
                ["ref", dict([["dist", "norm"], ["loc", 67], ["scale", 2]])],
                ["proposal", 2], ["latex", r"H_0"]])]])}],
    ["H_DESpriors", {
        _desc: "Hubble parameter (reduced range for DES and lensing-only constraints)",
        kinds.theory: {_camb: None, _classy: None},
        _params: dict([
            ["H0", dict([
                ["prior", dict([["min", 55], ["max", 91]])],
                ["ref", dict([["dist", "norm"], ["loc", 67], ["scale", 2]])],
                ["proposal", 2], ["latex", r"H_0"]])]])}],
    ["sound_horizon_last_scattering", {
        _desc: "Angular size of the sound horizon at last scattering "
               "(approximate, if using CAMB)",
        kinds.theory: {
            _camb: {
                _params: dict([
                    ["theta_MC_100", dict([
                        ["prior", dict([["min", 0.5], ["max", 10]])],
                        ["ref",
                         dict([["dist", "norm"], ["loc", 1.04109],
                               ["scale", 0.0004]])],
                        ["proposal", 0.0002],
                        ["latex", r"100\theta_\mathrm{MC}"],
                        ["drop", True], ["renames", "theta"]])],
                    ["cosmomc_theta", dict([
                        ["value", "lambda theta_MC_100: 1.e-2*theta_MC_100"],
                        ["derived", False]])],
                    ["H0", {"latex": r"H_0", "min": H0_min, "max": H0_max}]]),
                _extra_args: dict([["theta_H0_range", [H0_min, H0_max]]])},
            _classy: {
                _params: dict([
                    ["theta_s_1e2", dict([
                        ["prior", dict([["min", 0.5], ["max", 10]])],
                        ["ref", dict([
                            ["dist", "norm"], ["loc", 1.0416], ["scale", 0.0004]])],
                        ["proposal", 0.0002],
                        ["latex", r"100\theta_\mathrm{s}"],
                        ["drop", True]])],
                    ["100*theta_s", dict([
                        ["value", "lambda theta_s_1e2: theta_s_1e2"],
                        ["derived", False]])],
                    ["H0", {"latex": r"H_0"}]])}}}]])

# Matter sector (minus light species)
N_eff_std = 3.046
nu_mass_fac = 94.0708
matter: InfoDict = dict([
    ["omegab_h2, omegac_h2", {
        _desc: "Flat prior on Omega*h^2 for baryons and cold dark matter",
        kinds.theory: {_camb: None, _classy: None},
        _params: dict([
            ["omegabh2", dict([
                ["prior", dict([["min", 0.005], ["max", 0.1]])],
                ["ref",
                 dict([["dist", "norm"], ["loc", 0.0224], ["scale", 0.0001]])],
                ["proposal", 0.0001], ["latex", r"\Omega_\mathrm{b} h^2"]])],
            ["omegach2", dict([
                ["prior", dict([["min", 0.001], ["max", 0.99]])],
                ["ref",
                 dict([["dist", "norm"], ["loc", 0.120], ["scale", 0.001]])],
                ["proposal", 0.0005], ["latex", r"\Omega_\mathrm{c} h^2"]])],
            ["omegam", {"latex": r"\Omega_\mathrm{m}"}]])}],
    ["Omegab, Omegam", {
        _desc: "Flat prior on Omega for baryons and total matter",
        kinds.theory: {_camb: None, _classy: None},
        _params: dict([
            ["omegab", dict([
                ["prior", dict([["min", 0.03], ["max", 0.07]])],
                ["ref",
                 dict([["dist", "norm"], ["loc", 0.0495], ["scale", 0.004]])],
                ["proposal", 0.004], ["latex", r"\Omega_\mathrm{b}"],
                ["drop", True]])],
            ["omegam", dict([
                ["prior", dict([["min", 0.1], ["max", 0.9]])],
                ["ref",
                 dict([["dist", "norm"], ["loc", 0.316], ["scale", 0.02]])],
                ["proposal", 0.02], ["latex", r"\Omega_\mathrm{m}"],
                ["drop", True]])],
            ["omegabh2", dict([
                ["value", "lambda omegab, H0: omegab*(H0/100)**2"],
                ["latex", r"\Omega_\mathrm{b} h^2"]])],
            ["omegach2", dict([
                ["value", ("lambda omegam, omegab, mnu, H0: "
                           "(omegam-omegab)*(H0/100)**2-(mnu*(%g/3)**0.75)/%g" % (
                               N_eff_std, nu_mass_fac))],
                ["latex", r"\Omega_\mathrm{c} h^2"]])]
        ])}]])
for m in matter.values():
    m[_params]["omegamh2"] = dict([
        ["derived", "lambda omegam, H0: omegam*(H0/100)**2"],
        ["latex", r"\Omega_\mathrm{m} h^2"]])

# Neutrinos and other extra matter
neutrinos: InfoDict = dict([
    ["one_heavy_planck", {
        _desc: "Two massless nu and one with m=0.06. Neff=3.046",
        kinds.theory: {
            _camb: {_extra_args: {"num_massive_neutrinos": 1, "nnu": 3.046},
                    _params: dict([["mnu", 0.06]])},
            _classy: {_extra_args: {"N_ncdm": 1, "N_ur": 2.0328},
                      _params: dict([
                          ["m_ncdm", {"value": 0.06, "renames": "mnu"}]])}}}],
    ["varying_mnu", {
        _desc: "Varying total mass of 3 degenerate nu's, with N_eff=3.046",
        kinds.theory: {
            _camb: {_extra_args: {"num_massive_neutrinos": 3, "nnu": 3.046},
                    _params: dict([
                        ["mnu", dict([
                            ["prior", dict([["min", 0], ["max", 5]])],
                            ["ref", dict([
                                ["dist", "norm"], ["loc", 0.02], ["scale", 0.1]])],
                            ["proposal", 0.03], ["latex", r"\sum m_\nu"]])]])},
            _classy: {_extra_args: {"N_ncdm": 1, "deg_ncdm": 3, "N_ur": 0.00641},
                      _params: dict([
                          ["m_ncdm", dict([
                              ["prior", dict([["min", 0], ["max", 1.667]])],
                              ["ref", dict([
                                  ["dist", "norm"], ["loc", 0.0067],
                                  ["scale", 0.033]])],
                              ["proposal", 0.01], ["latex", r"m_\nu"]])],
                          ["mnu", dict([["derived", "lambda m_ncdm: 3 * m_ncdm"],
                                        ["latex", r"\sum m_\nu"]])]])}}}],
    ["varying_Neff", {
        _desc: "Varying Neff with two massless nu and one with m=0.06",
        kinds.theory: {
            _camb: {_extra_args: {"num_massive_neutrinos": 1},
                    _params: dict([
                        ["mnu", 0.06],
                        ["nnu", dict([
                            ["prior", dict([["min", 0.05], ["max", 10]])],
                            ["ref", dict([
                                ["dist", "norm"], ["loc", 3.046], ["scale", 0.05]])],
                            ["proposal", 0.05],
                            ["latex", r"N_\mathrm{eff}"]])]])},
            _classy: {_extra_args: {"N_ncdm": 1},
                      _params: dict([
                          ["m_ncdm",
                           dict([["value", 0.06], ["renames", "mnu"]])],
                          ["N_ur", dict([
                              ["prior", dict([["min", 0.0001], ["max", 9]])],
                              ["ref", dict([
                                  ["dist", "norm"], ["loc", 2.0328],
                                  ["scale", 0.05]])],
                              ["proposal", 0.05],
                              ["latex", r"N_\mathrm{ur}"]])],
                          ["nnu", dict([
                              ["derived", "lambda Neff: Neff"],
                              ["latex", r"N_\mathrm{eff}"]])]])}}}]])
neutrinos.update(dict([
    ["varying_mnu_Neff", {
        _desc: "Varying Neff and total mass of 3 degenerate nu's",
        kinds.theory: {
            _camb: {
                _extra_args: {"num_massive_neutrinos": 3},
                _params: dict([
                    ["mnu", deepcopy(
                        neutrinos["varying_mnu"][kinds.theory]["camb"][_params]["mnu"])],
                    ["nnu", deepcopy(
                        neutrinos["varying_Neff"][kinds.theory]["camb"][_params][
                            "nnu"])]])},
            #            _classy: {
            #                _extra_args: {"N_ncdm": 1, "deg_ncdm": 3},
            #                _params: dict([
            #                    ["m_ncdm", deepcopy(
            #                        neutrinos["varying_mnu"][_theory]["classy"][_params]["m_ncdm"])],
            #                    ["mnu", deepcopy(
            #                        neutrinos["varying_mnu"][_theory]["classy"][_params]["mnu"])],
            #                    ["N_ur", deepcopy(
            #                        neutrinos["varying_Neff"][_theory]["classy"][_params]["N_ur"])],
            #                    ["nnu", deepcopy(
            #                        neutrinos["varying_Neff"][_theory]["classy"][_params]["nnu"])]
            #                ])}
        }}]]))
# neutrinos["varying_mnu_Neff"][_theory][_classy][_params]["N_ur"]["ref"]["loc"] = 0.00641

#    # ["varying_Neff+1sterile", {
#    #     _desc: "Varying Neff plus 1 sterile neutrino (SM nu's with m=0,0,0.06)",
#    #     _theory: {
#    #         _camb: {_extra_args: dict(
#    #             list(neutrinos["varying_Neff"][_theory][_camb][_extra_args].items()) +
#    #             [["accuracy_level", 1.2]])}},
#    #     _params: dict([
#    #         ["nnu", dict([
#    #             ["prior", dict([["min", 3.046], ["max", 10]])],
#    #             ["ref", dict([["dist", "norm"], ["loc", 3.046], ["scale", 0.05]])],
#    #             ["proposal", 0.05], ["latex", r"N_\mathrm{eff}"]])],
#    #         ["meffsterile", dict([
#    #             ["prior", dict([["min", 0], ["max", 3]])#,
#    #             ["ref", dict([["dist", "norm"], ["loc", 0.1], ["scale", 0.1]])],
#    #             ["proposal", 0.03],
#    #             ["latex", r"m_{\nu,\mathrm{sterile}}^\mathrm{eff}"]])]])}]

# Dark Energy
dark_energy: InfoDict = dict([
    ["lambda", {
        _desc: "Cosmological constant (w=-1)",
        kinds.theory: {_camb: None, _classy: None},
        _params: dict([
            ["omegal", {"latex": r"\Omega_\Lambda"}]])}],
    ["de_w", {
        _desc: "Varying constant eq of state",
        kinds.theory: {_camb: None,
                       _classy: {_params: {"Omega_Lambda": 0}}},
        _params: dict([
            ["w", dict([
                ["prior", dict([["min", -3], ["max", -0.333]])],
                ["ref",
                 dict([["dist", "norm"], ["loc", -0.99], ["scale", 0.02]])],
                ["proposal", 0.02], ["latex", r"w_\mathrm{DE}"]])]])}],
    ["de_w_wa", {
        _desc: "Varying constant eq of state with w(a) = w0 + (1-a) wa",
        kinds.theory: {_camb: {_extra_args: {'dark_energy_model': 'ppf'}},
                       _classy: {_params: {"Omega_Lambda": 0}}},
        _params: dict([
            ["w", dict([
                ["prior", dict([["min", -3], ["max", 1]])],
                ["ref",
                 dict([["dist", "norm"], ["loc", -0.99], ["scale", 0.02]])],
                ["proposal", 0.02], ["latex", r"w_{0,\mathrm{DE}}"]])],
            ["wa", dict([
                ["prior", dict([["min", -3], ["max", 2]])],
                ["ref", dict([["dist", "norm"], ["loc", 0], ["scale", 0.05]])],
                ["proposal", 0.05], ["latex", r"w_{a,\mathrm{DE}}"]])]])}]])

# BBN
bbn_derived_camb: InfoDict = dict([
    ["YpBBN", dict([["latex", r"Y_P^\mathrm{BBN}"]])],
    ["DHBBN", dict([["derived", r"lambda DH: 10**5*DH"],
                    ["latex", r"10^5 \mathrm{D}/\mathrm{H}"]])]])
bbn = dict([
    ["consistency", {
        _desc: "Primordial Helium fraction inferred from BBN consistency",
        kinds.theory: {_camb: {_params: bbn_derived_camb},
                       _classy: None},
        _params: dict([
            ["yheused", {"latex": r"Y_\mathrm{P}"}]])}],
    ["YHe_des_y1", {
        _desc: "Fixed Y_P = 0.245341 (used in DES Y1)",
        kinds.theory: {_camb: None,
                       _classy: None},
        _params: dict([
            ["yhe", 0.245341]])}],
    ["YHe", {
        _desc: "Varying primordial Helium fraction",
        kinds.theory: {_camb: None,
                       _classy: None},
        _params: dict([
            ["yhe", dict([
                ["prior", dict([["min", 0.1], ["max", 0.5]])],
                ["ref",
                 dict([["dist", "norm"], ["loc", 0.237], ["scale", 0.006]])],
                ["proposal", 0.006], ["latex", r"Y_\mathrm{P}"]])]])}], ])

# Reionization
reionization: InfoDict = dict([
    ["std", {
        _desc: "Standard reio, lasting delta_z=0.5",
        kinds.theory: {_camb: None, _classy: None},
        _params: dict([
            ["tau", dict([
                ["prior", dict([["min", 0.01], ["max", 0.8]])],
                ["ref",
                 dict([["dist", "norm"], ["loc", 0.055], ["scale", 0.006]])],
                ["proposal", 0.003], ["latex", r"\tau_\mathrm{reio}"]])],
            ["zrei", {"latex": r"z_\mathrm{re}"}]])}],
    ["gauss_prior", {
        _desc: "Standard reio, lasting delta_z=0.5, gaussian prior around tau=0.07",
        kinds.theory: {_camb: None, _classy: None},
        _params: dict([
            ["tau", dict([
                ["prior",
                 dict([["dist", "norm"], ["loc", 0.07], ["scale", 0.02]])],
                ["ref",
                 dict([["dist", "norm"], ["loc", 0.07], ["scale", 0.01]])],
                ["proposal", 0.005], ["latex", r"\tau_\mathrm{reio}"]])],
            ["zrei", {"latex": r"z_\mathrm{re}"}]])}],
    ["irrelevant", {
        _desc: "Irrelevant (NB: only valid for non-CMB or CMB-marged datasets!)",
        kinds.theory: {_camb: None, _classy: None},
        _params: dict([])}], ])

# EXPERIMENTS ############################################################################
base_precision: InfoDict = {_camb: {"halofit_version": "mead"},
                            _classy: {"non linear": "hmcode", "hmcode_min_k_max": 20}}
cmb_precision = deepcopy(base_precision)
cmb_precision[_camb].update({"bbn_predictor": "PArthENoPE_880.2_standard.dat",
                             "lens_potential_accuracy": 1})
cmb_sampler_recommended: InfoDict = {"mcmc": {
    "drag": True, "oversample_power": 0.4, "proposal_scale": 1.9}}

like_cmb: InfoDict = dict([
    [_none, {}],
    ["planck_2018", {
        _desc: "Planck 2018 (Polarized CMB + lensing)",
        _comment: None,
        kinds.sampler: cmb_sampler_recommended,
        kinds.theory: {theo: {_extra_args: cmb_precision[theo]}
                       for theo in [_camb, _classy]},
        kinds.likelihood: dict([
            ["planck_2018_lowl.TT", None],
            ["planck_2018_lowl.EE", None],
            ["planck_2018_highl_plik.TTTEEE", None],
            ["planck_2018_lensing.clik", None]])}],
    ["planck_2018_bk15", {
        _desc: "Planck 2018 (Polarized CMB + lensing) + Bicep/Keck-Array 2015",
        kinds.sampler: cmb_sampler_recommended,
        kinds.theory: {theo: {_extra_args: cmb_precision[theo]}
                       for theo in [_camb, _classy]},
        kinds.likelihood: dict([
            ["planck_2018_lowl.TT", None],
            ["planck_2018_lowl.EE", None],
            ["planck_2018_highl_plik.TTTEEE", None],
            ["planck_2018_lensing.clik", None],
            ["bicep_keck_2015", None]])}],
    ["planck_2018_CMBmarged_lensing", {
        _desc: "Planck 2018 CMB-marginalized lensing",
        kinds.sampler: cmb_sampler_recommended,
        kinds.theory: {theo: {_extra_args: cmb_precision[theo]}
                       for theo in [_camb, _classy]},
        kinds.likelihood: dict([
            ["planck_2018_lensing.CMBMarged", None]])}],
])
like_cmb["planck_2018_bk15"][_comment] = like_cmb["planck_2018"][_comment]
# Add common CMB derived parameters
for name, m in like_cmb.items():
    # Don't add the derived parameter to the no-CMB case!
    if not m:
        continue
    if _params not in m:
        m[_params] = dict()
    m[_params].update(dict([
        ["sigma8", {"latex": r"\sigma_8"}],
        ["s8h5", dict([
            ["derived", "lambda sigma8, H0: sigma8*(H0*1e-2)**(-0.5)"],
            ["latex", r"\sigma_8/h^{0.5}"]])],
        ["s8omegamp5", dict([
            ["derived", "lambda sigma8, omegam: sigma8*omegam**0.5"],
            ["latex", r"\sigma_8 \Omega_\mathrm{m}^{0.5}"]])],
        ["s8omegamp25", dict([
            ["derived", "lambda sigma8, omegam: sigma8*omegam**0.25"],
            ["latex", r"\sigma_8 \Omega_\mathrm{m}^{0.25}"]])],
        ["A", dict([
            ["derived", "lambda As: 1e9*As"],
            ["latex", r"10^9 A_\mathrm{s}"]])],
        ["clamp", dict([
            ["derived", "lambda As, tau: 1e9*As*np.exp(-2*tau)"],
            ["latex", r"10^9 A_\mathrm{s} e^{-2\tau}"]])],
        ["age", dict([
            ["latex", r"{\rm{Age}}/\mathrm{Gyr}"]])],
        ["rdrag", dict([
            ["latex", r"r_\mathrm{drag}"]])]]))
    if "cmbmarged" in name.lower():
        m[_params].pop("A")
        m[_params].pop("clamp")
# Some more, in case we want to add them at some point, described in
# https://wiki.cosmos.esa.int/planckpla2015/images/b/b9/Parameter_tag_definitions_2015.pdf
#    "zstar":       {"latex": r"z_*"},
#    "rstar":       {"latex": r"r_*"},
#    "thetastar":   {"latex": r"100\theta_*"},
#    "DAstar":      {"latex": r"D_\mathrm{A}/\mathrm{Gpc}"},
#    "zdrag":       {"latex": r"z_\mathrm{drag}"},
#    "kd":          {"latex": r"k_\mathrm{D}"},
#    "thetad":      {"latex": r"100\theta_\mathrm{D}"},
#    "zeq":         {"latex": r"z_\mathrm{eq}"},
#    "keq":         {"latex": r"k_\mathrm{eq}"},
#    "thetaeq":     {"latex": r"100\theta_\mathrm{eq}"},
#    "thetarseq":   {"latex": r"100\theta_\mathrm{s,eq}"},

like_bao: InfoDict = dict([
    [_none, {}],
    ["BAO_planck_2018", {
        _desc: "Baryon acoustic oscillation data from DR12, MGS and 6DF",
        kinds.theory: {_camb: None, _classy: None},
        kinds.likelihood: dict([
            ["bao.sixdf_2011_bao", None],
            ["bao.sdss_dr7_mgs", None],
            ["bao.sdss_dr12_consensus_bao", None]])}],
])

like_des: InfoDict = dict([
    [_none, {}],
    ["des_y1_clustering", {
        _desc: "Galaxy clustering from DES Y1",
        kinds.theory: {theo: {_extra_args: base_precision[theo]}
                       for theo in [_camb, _classy]},
        kinds.likelihood: dict([
            ["des_y1.clustering", None]])}],
    ["des_y1_galaxy_galaxy", {
        _desc: "Galaxy-galaxy lensing from DES Y1",
        kinds.theory: {theo: {_extra_args: base_precision[theo]}
                       for theo in [_camb, _classy]},
        kinds.likelihood: dict([
            ["des_y1.galaxy_galaxy", None]])}],
    ["des_y1_shear", {
        _desc: "Cosmic shear data from DES Y1",
        kinds.theory: {theo: {_extra_args: base_precision[theo]}
                       for theo in [_camb, _classy]},
        kinds.likelihood: dict([
            ["des_y1.shear", None]])}],
    ["des_y1_joint", {
        _desc: "Combination of galaxy clustering and weak lensing data from DES Y1",
        kinds.theory: {theo: {_extra_args: base_precision[theo]}
                       for theo in [_camb, _classy]},
        kinds.likelihood: dict([
            ["des_y1.joint", None]])}],
])

like_sn: InfoDict = dict([
    [_none, {}],
    ["Pantheon", {
        _desc: "Supernovae data from the Pantheon sample",
        kinds.theory: {_camb: None, _classy: None},
        kinds.likelihood: dict([
            ["sn.pantheon", None]])}],
])

like_H0: InfoDict = dict([
    [_none, {}],
    ["Riess2018a", {
        _desc: "Local H0 measurement from Riess et al. 2018a (used in Planck 2018)",
        kinds.theory: {_camb: None, _classy: None},
        kinds.likelihood: dict([
            ["H0.riess2018a", None]])}],
    ["Riess201903", {
        _desc: "Local H0 measurement from Riess et al. 2019",
        kinds.theory: {_camb: None, _classy: None},
        kinds.likelihood: dict([
            ["H0.riess201903", None]])}],
])

# SAMPLERS ###############################################################################

sampler: InfoDict = dict([
    ["MCMC", {
        _desc: "MCMC sampler with covmat learning",
        kinds.sampler: {"mcmc": {"covmat": "auto"}}}],
    ["PolyChord", {
        _desc: "Nested sampler, affine invariant and multi-modal",
        kinds.sampler: {"polychord": None}}], ])

# PRESETS ################################################################################

planck_base_model = {
    "primordial": "SFSR",
    "geometry": "flat",
    "hubble": "sound_horizon_last_scattering",
    "matter": "omegab_h2, omegac_h2",
    "neutrinos": "one_heavy_planck",
    "dark_energy": "lambda",
    "bbn": "consistency",
    "reionization": "std"}
default_sampler = {"sampler": "MCMC"}

preset: InfoDict = dict([
    [_none, {_desc: "(No preset chosen)"}],
    # Pure CMB #######################################################
    ["planck_2018_camb", {
        _desc: "Planck 2018 with CAMB",
        "theory": _camb,
        "like_cmb": "planck_2018"}],
    ["planck_2018_classy", {
        _desc: "Planck 2018 with CLASS",
        "theory": _classy,
        "like_cmb": "planck_2018"}],
    ["planck_2018_bicep_camb", {
        _desc: "Planck 2018 + BK15 (with tensor modes) with CAMB",
        "theory": _camb,
        "primordial": "SFSR_t",
        "like_cmb": "planck_2018_bk15"}],
    ["planck_2018_bicep_classy", {
        _desc: "Planck 2018 + BK15 (with tensor modes) with CLASS",
        "theory": _classy,
        "primordial": "SFSR_t",
        "like_cmb": "planck_2018_bk15"}],
    # CMB+BAO ######################################################
    ["planck_2018_BAO_camb", {
        _desc: "Planck 2018 + BAO with CAMB",
        "theory": _camb,
        "like_cmb": "planck_2018",
        "like_bao": "BAO_planck_2018"}],
    ["planck_2018_BAO_classy", {
        _desc: "Planck 2018 + BAO with CLASS",
        "theory": _classy,
        "like_cmb": "planck_2018",
        "like_bao": "BAO_planck_2018"}],
    # CMB+BAO+SN ###################################################
    ["planck_2018_BAO_SN_camb", {
        _desc: "Planck 2018 + BAO + SN with CAMB",
        "theory": _camb,
        "like_cmb": "planck_2018",
        "like_bao": "BAO_planck_2018",
        "like_sn": "Pantheon"}],
    ["planck_2018_BAO_SN_classy", {
        _desc: "Planck 2018 + BAO + SN with CLASS",
        "theory": _classy,
        "like_cmb": "planck_2018",
        "like_bao": "BAO_planck_2018",
        "like_sn": "Pantheon"}],
    # CMB+DES+BAO+SN ###################################################
    ["planck_2018_DES_BAO_SN_camb", {
        _desc: "Planck 2018 + DESjoint + BAO + SN with CAMB",
        "theory": _camb,
        "like_cmb": "planck_2018",
        "like_bao": "BAO_planck_2018",
        "like_des": "des_y1_joint",
        "like_sn": "Pantheon"}],
    ["planck_2018_DES_BAO_SN_classy", {
        _desc: "Planck 2018 + DESjoint + BAO + SN with CLASS",
        "theory": _classy,
        "like_cmb": "planck_2018",
        "like_bao": "BAO_planck_2018",
        "like_des": "des_y1_joint",
        "like_sn": "Pantheon"}],
])

# Add planck baseline model
for pre in preset.values():
    pre.update(
        {field: value for field, value in planck_base_model.items() if field not in pre})
    pre.update(default_sampler)

# Lensing-only ###################################################
preset.update(dict([
    [_none, {_desc: "(No preset chosen)"}],
    ["planck_2018_DES_lensingonly_camb", {
        _desc: "Planck 2018 + DES Y1 lensing-only with CAMB",
        "theory": _camb,
        "like_cmb": "planck_2018_CMBmarged_lensing",
        "like_des": "des_y1_shear"}],
    ["planck_2018_DES_lensingonly_classy", {
        _desc: "Planck 2018 + DES Y1 lensing-only with CLASS",
        "theory": _classy,
        "like_cmb": "planck_2018_CMBmarged_lensing",
        "like_des": "des_y1_shear"}],
]))

lensingonly_model = {
    "primordial": "SFSR_DESpriors",
    "geometry": "flat",
    "hubble": "H_DESpriors",
    "matter": "Omegab, Omegam",
    "neutrinos": "one_heavy_planck",
    "dark_energy": "lambda",
    "bbn": "YHe_des_y1",
    "reionization": "irrelevant"}

# Add planck baseline model
for name, pre in preset.items():
    if "lensingonly" in name:
        pre.update(
            {field: value for field, value in lensingonly_model.items() if
             field not in pre})
        pre.update(default_sampler)

# BASIC INSTALLATION #####################################################################
install_basic: InfoDict = {
    kinds.theory: {_camb: None, _classy: None},
    kinds.likelihood: {
        "planck_2018_lowl.TT": None,
        "planck_2018_lensing.native": None,
        "bicep_keck_2015": None,
        "sn.pantheon": None,
        "bao.sdss_dr12_consensus_final": None,
        "des_y1.joint": None}}

install_tests = deepcopy(install_basic)
install_tests[kinds.likelihood].update({"planck_2015_lowl": None,
                                        "planck_2018_highl_plik.TT_unbinned": None,
                                        "planck_2018_highl_plik.TT_lite_native": None,
                                        "planck_2018_highl_CamSpec.TT": None,
                                        "planck_2018_highl_CamSpec.TT_native": None,
                                        })

# CONTENTS FOR COMBO-BOXED IN A GUI ######################################################

_combo_dict_text = (
    ["Presets", (["preset", "Presets"],)],
    ["Cosmological Model", (
        ["theory", "Theory code"],
        ["primordial", "Primordial perturbations"],
        ["geometry", "Geometry"],
        ["hubble", "Hubble parameter constraint"],
        ["matter", "Matter sector"],
        ["neutrinos", "Neutrinos and other extra matter"],
        ["dark_energy", "Lambda / Dark energy"],
        ["bbn", "BBN"],
        ["reionization", "Reionization history"])],
    ["Data sets", (
        ["like_cmb", "CMB experiments"],
        ["like_bao", "BAO experiments"],
        ["like_des", "DES measurements"],
        ["like_sn", "SN experiments"],
        ["like_H0", "Local H0 measurements"])],
    ["Sampler", (["sampler", "Samplers"],)])
