import os
from copy import deepcopy
from flaky import flaky
from scipy.stats import multivariate_normal
from getdist.mcsamples import loadMCSamples, MCSamplesFromCobaya
import numpy as np
import pytest

from cobaya.run import run
from cobaya.post import post
from cobaya.tools import KL_norm
from cobaya.typing import ParamsDict, InputDict
from cobaya.conventions import _output_prefix, _params, _force, kinds
from cobaya.conventions import _prior, partag, _separator_files
from cobaya.conventions import _post, _post_add, _post_remove, _post_suffix
from cobaya import mpi

pytestmark = pytest.mark.mpi

_post_ = _separator_files + _post + _separator_files

mean = np.array([0, 0])
sigma = 0.5
cov = np.array([[sigma ** 2, 0], [0, sigma ** 2]])
sampled = {"mean": mean, "cov": cov}
target = {"mean": mean + np.array([sigma / 2, 0]), "cov": cov}


def sampled_pdf(a, b):
    return multivariate_normal.logpdf([a, b], mean=sampled["mean"], cov=sampled["cov"])


def target_pdf(a, b, c=0):
    logp = multivariate_normal.logpdf([a, b], mean=target["mean"], cov=target["cov"])
    derived = {"cprime": c}
    return logp, derived


target_pdf_prior = lambda a, b, c=0: target_pdf(a, b, c=0)[0]

_range = {"min": -2, "max": 2}
ref_pdf = {"dist": "norm", "loc": 0, "scale": 0.1}
info_params: ParamsDict = dict([
    ("a", {"prior": _range, "ref": ref_pdf, partag.proposal: sigma}),
    ("b", {"prior": _range, "ref": ref_pdf, partag.proposal: sigma}),
    ("a_plus_b", {partag.derived: lambda a, b: a + b})])

info_sampler = {"mcmc": {"Rminus1_stop": 0.005}}
info_sampler_dummy = {"evaluate": {"N": 10}}


@flaky(max_runs=3, min_passes=1)
@mpi.sync_errors
def test_post_prior(tmpdir):
    # Generate original chain
    info: InputDict = {
        _output_prefix: os.path.join(tmpdir, "gaussian"), _force: True,
        _params: info_params, kinds.sampler: info_sampler,
        kinds.likelihood: {"one": None}, _prior: {"gaussian": sampled_pdf}}
    info_post: InputDict = {
        _output_prefix: info[_output_prefix], _force: True,
        _post: {_post_suffix: "foo", 'skip': 0.05,
                _post_remove: {_prior: {"gaussian": None}},
                _post_add: {_prior: {"target": target_pdf_prior}}}}
    _, sampler = run(info)
    for mem in [False, True]:
        post(info_post, sample=sampler.products()["sample"] if mem else None)

        # Load with GetDist and compare
        if mpi.is_main_process():
            mcsamples = loadMCSamples(
                info_post[_output_prefix] + _post_ + info_post[_post][_post_suffix])
            new_mean = mcsamples.mean(["a", "b"])
            new_cov = mcsamples.getCovMat().matrix
            mpi.share((new_mean, new_cov))
        else:
            new_mean, new_cov = mpi.share()
        assert abs(KL_norm(target["mean"], target["cov"], new_mean, new_cov)) < 0.02


@flaky(max_runs=3, min_passes=1)
def test_post_likelihood():
    """
    Swaps likelihood "gaussian" for "target".

    It also tests aggregated chi2's by removing and adding a likelihood to an existing
    type.
    """
    # Generate original chain
    info_params_local = deepcopy(info_params)
    info_params_local["dummy"] = 0
    dummy_loglike_add = 0.1
    dummy_loglike_remove = 0.01
    info = {
        _output_prefix: None, _force: True,
        _params: info_params_local, kinds.sampler: info_sampler,
        kinds.likelihood: {
            "gaussian": {"external": sampled_pdf, "type": "AA"},
            "dummy": {"external": lambda dummy: 1, "type": "BB"},
            "dummy_remove": {"external": lambda dummy: dummy_loglike_add, "type": "BB"}}}
    info_out, sampler = run(info)
    info_out.update({
        _post: {_post_suffix: "foo",
                _post_remove: {kinds.likelihood: {
                    "gaussian": None, "dummy_remove": None}},
                _post_add: {kinds.likelihood: {
                    "target": {
                        "external": target_pdf, "type": "AA",
                        "output_params": ["cprime"]},
                    "dummy_add": {
                        "external": lambda dummy: dummy_loglike_remove, "type": "BB"}}}}})
    info_post_out, products_post = post(info_out, sampler.products()["sample"])
    samples = mpi.gather(products_post["sample"])

    # Load with GetDist and compare
    if mpi.is_main_process():
        mcsamples = MCSamplesFromCobaya(info_post_out, samples, name_tag="sample")
        new_mean = mcsamples.mean(["a", "b"])
        new_cov = mcsamples.getCovMat().matrix
        mpi.share((new_mean, new_cov))
    else:
        new_mean, new_cov = mpi.share()
    assert abs(KL_norm(target["mean"], target["cov"], new_mean, new_cov)) < 0.02
    assert np.allclose(products_post["sample"]["chi2__AA"],
                       products_post["sample"]["chi2__target"])
    assert np.allclose(products_post["sample"]["chi2__BB"],
                       products_post["sample"]["chi2__dummy"] +
                       products_post["sample"]["chi2__dummy_add"])


def test_post_params():
    # Tests:
    # - added simple dynamical derived parameter "a+b"
    # - added dynamical derived parameter that depends on *new* chi2__target
    # - added new fixed input "c" + new derived-from-external-function "cprime"
    # Generate original chain
    info = {
        _params: info_params, kinds.sampler: info_sampler_dummy,
        kinds.likelihood: {"gaussian": sampled_pdf}}
    updated_info_gaussian, sampler_gaussian = run(info)
    products_gaussian = sampler_gaussian.products()
    info_post = {
        _post: {_post_suffix: "foo",
                _post_remove: {_params: {"a_plus_b": None}},
                _post_add: {
                    kinds.likelihood: {
                        "target": {"external": target_pdf, "output_params": ["cprime"]}},
                    _params: {
                        "c": 1.234,
                        "a_minus_b": {partag.derived: "lambda a,b: a-b"},
                        "my_chi2__target": {
                            partag.derived: "lambda chi2__target: chi2__target"},
                        "cprime": None}}}}
    info_post.update(updated_info_gaussian)
    products = post(info_post, products_gaussian["sample"]).products
    # Compare parameters
    assert np.allclose(
        products["sample"]["a"] - products["sample"]["b"],
        products["sample"]["a_minus_b"])
    assert np.allclose(
        products["sample"]["cprime"], info_post[_post][_post_add][_params]["c"])
    assert np.allclose(
        products["sample"]["my_chi2__target"], products["sample"]["chi2__target"])
