#!usr/bin/python
import copy
from .constants import Model, FungicideEffectType
from functools import reduce
from .fungicide import Fungicide, FungicideEffect
import itertools
from . import types
import typing

BetaType = typing.Union[types.Number, typing.List[types.Number]]


class BetaMap(object):
    def __init__(self) -> None:
        self.baseStore = {}  # type: typing.Any

        for i, j, k, l in itertools.product((True, False), repeat = 4):
            lastDict = self.baseStore
            for m in (i, j, k):
                if m not in lastDict:
                    lastDict[m] = {}
                lastDict = lastDict[m]
            self.baseStore[i][j][k][l] = None

        self.store = {}  # type: typing.Any

    def addModel(self, model: str) -> None:
        self.store[model] = copy.deepcopy(self.baseStore)

    def get(self, model: str, densityDependentInfection: bool, seasonality: bool, latentPeriod: bool, fungicideDecay: bool) -> BetaType:
        return self.store[model][densityDependentInfection][seasonality][latentPeriod][fungicideDecay]

    def set(self, model: str, densityDependentInfection: bool, seasonality: bool, latentPeriod: bool, fungicideDecay: bool, value: BetaType) -> None:
        self.store[model][densityDependentInfection][seasonality][latentPeriod][fungicideDecay] = value


class Parameters(object):
    def __init__(self, model: str) -> None:
        if model == Model.WheatSeptoria:
            self.setToWheatSeptoria()
        elif model == Model.GrapePowderyMildew:
            self.setToGrapePowderyMildew()
        else:
            raise TypeError

    @property
    def densityDependentInfection(self) -> bool:
        return self._densityDependentInfection

    @densityDependentInfection.setter
    def densityDependentInfection(self, value: bool) -> None:
        self._densityDependentInfection = value
        self._setFittedParams()

    @property
    def fungicideDecay(self) -> bool:
        return self._fungicideDecay

    @fungicideDecay.setter
    def fungicideDecay(self, value: bool) -> None:
        self._fungicideDecay = value
        self._setFittedParams()

    @property
    def seasonality(self) -> bool:
        return self._seasonality

    @seasonality.setter
    def seasonality(self, value: bool) -> None:
        self._seasonality = value
        self._setFittedParams()

    @property
    def latentPeriod(self) -> bool:
        return self._latentPeriod

    @latentPeriod.setter
    def latentPeriod(self, value: bool) -> None:
        self._latentPeriod = value
        self._setFittedParams()

    @property
    def maxSeasons(self) -> types.Number:
        if self.model == Model.WheatSeptoria:
            val = self.maxTime / (self.GS87 - self.seasonStartTime)
        elif self.model == Model.GrapePowderyMildew:
            val = self.maxTime / (self.seasonEndTime - self.budBreakTime)

        return val

    @maxSeasons.setter
    def maxSeasons(self, value: types.Number) -> None:
        if self.model == Model.WheatSeptoria:
            self.maxTime = value * (self.GS87 - self.seasonStartTime)
        elif self.model == Model.GrapePowderyMildew:
            self.maxTime = value * (self.seasonEndTime - self.budBreakTime)

    def _setFittedParams(self) -> None:
        self.beta = self.betaMap.get(
            self.model,
            self.densityDependentInfection,
            self.seasonality,
            self.latentPeriod,
            self.fungicideDecay
        )

        if self.model == Model.GrapePowderyMildew:
            if not self.seasonality:
                self.r = [0.0923]  # From fit
                self.k = [94810]  # From fit
            else:
                self.r = [0.147, 0.032]  # Burie paper
                self.k = [26106, 2.461E8]  # Pers. comm Burie

    betaMap = BetaMap()

    @property
    def ER0(self) -> types.Number:
        if self.latentPeriod:
            if self.model == Model.GrapePowderyMildew:
                return self.phi * self.rFreq
            elif self.model == Model.WheatSeptoria and not self.seasonality:
                return self.phi * self.rFreq
            else:
                return self._ER0
        else:
            return 0

    @ER0.setter
    def ER0(self, value: types.Number) -> None:
        if self.latentPeriod:
            if self.model == Model.GrapePowderyMildew:
                self.phi = value / self.rFreq
            elif self.model == Model.WheatSeptoria and not self.seasonality:
                self.phi = value / self.rFreq
            else:
                self._ER0 = value

    @property
    def ES0(self) -> types.Number:
        if self.latentPeriod:
            if self.model == Model.GrapePowderyMildew:
                return self.phi * (1 - self.rFreq)
            elif self.model == Model.WheatSeptoria and not self.seasonality:
                return self.phi * (1 - self.rFreq)
            else:
                return self._ES0
        else:
            return 0

    @ES0.setter
    def ES0(self, value: types.Number) -> None:
        if self.latentPeriod:
            if self.model == Model.GrapePowderyMildew:
                self.phi = value / (1 - self.rFreq)
            elif self.model == Model.WheatSeptoria and not self.seasonality:
                self.phi = value / (1 - self.rFreq)
            else:
                self._ES0 = value

    @property
    def PR0(self) -> types.Number:
        if self.model == Model.WheatSeptoria and self.seasonality:
            return self.rFreq * self.phi
        else:
            return 0

    @PR0.setter
    def PR0(self, value: types.Number) -> None:
        if self.model == Model.WheatSeptoria:
            self.phi = value / self.rFreq

    @property
    def PS0(self) -> types.Number:
        if self.model == Model.WheatSeptoria and self.seasonality:
            return self.phi * (1 - self.rFreq)
        else:
            return 0

    @PS0.setter
    def PS0(self, value: types.Number) -> None:
        if self.model == Model.WheatSeptoria:
            self.phi = value / (1 - self.rFreq)

    @property
    def IR0(self) -> types.Number:
        if self.model == Model.WheatSeptoria:
            if not self.seasonality and not self.latentPeriod:
                return self.rFreq * self.phi
            else:
                return self._IR0
        elif self.model == Model.GrapePowderyMildew and not self.latentPeriod:
                return self.rFreq * self.phi
        return self._IR0

    @IR0.setter
    def IR0(self, value: types.Number) -> None:
        if self.model == Model.WheatSeptoria:
            if not self.seasonality and not self.latentPeriod:
                self.phi = value / self.rFreq
            else:
                self._IR0 = value
        elif self.model == Model.GrapePowderyMildew and not self.latentPeriod:
                self.phi = value / self.rFreq
        self._IR0 = value

    @property
    def IS0(self) -> types.Number:
        if self.model == Model.WheatSeptoria:
            if not self.seasonality and not self.latentPeriod:
                return (1 - self.rFreq) * self.phi
            else:
                return self._IR0
        elif self.model == Model.GrapePowderyMildew and not self.latentPeriod:
                return (1 - self.rFreq) * self.phi
        return self._IR0

    @IS0.setter
    def IS0(self, value: types.Number) -> None:
        if self.model == Model.WheatSeptoria:
            if not self.seasonality and not self.latentPeriod:
                self.phi = value / (1 - self.rFreq)
            else:
                self._IS0 = value
        elif self.model == Model.GrapePowderyMildew and not self.latentPeriod:
                self.phi = value / (1 - self.rFreq)
        self._IS0 = value

    def setToWheatSeptoria(self) -> None:
        self._clearAttributes()

        self.model = Model.WheatSeptoria
        self.gamma = 1 / 266
        self.mu = 1 / 456
        self.highRisk = Fungicide(
            1.11E-2,
            **{
                FungicideEffectType.Protectant: FungicideEffect(1, 9.6),
                FungicideEffectType.Eradicant: FungicideEffect(1, 9.6)
            }
        )
        self.lowRisk = Fungicide(
            6.91E-3,
            **{
                FungicideEffectType.Protectant: FungicideEffect(0.48, 9.9)
            }
        )
        self.psi = 1E-10
        self.phi = 0.011
        self.S0 = 0.05
        self._ER0 = 0
        self._ES0 = 0
        self._IS0 = 0
        self._IR0 = 0
        self.R0 = 0

        # From Femke
        self.k = 4.2

        # From Peter
        self.r = 1.26E-2

        self.nu = 8.5E-3
        self.GS32 = 1456
        self.GS39 = 1700
        self.GS61 = 2066
        self.GS87 = 2900
        self.criticalYieldLoss = 0.05
        self.seasonStartTime = self.GS32 - 2 * 122  # Emergence of leaf 5
        self.maxTime = self.GS87 - self.seasonStartTime
        self.rFreq = self.psi
        self._densityDependentInfection = True
        self._fungicideDecay = True
        self._seasonality = True
        self._latentPeriod = True

        self._setFittedParams()

    def setToGrapePowderyMildew(self) -> None:
        self._clearAttributes()

        self.model = Model.GrapePowderyMildew
        self.gamma = 1 / 10
        self.mu = 1 / 10
        self.highRisk = Fungicide(
            0.231,
            **{
                FungicideEffectType.Protectant: FungicideEffect(1, 12.4),
                FungicideEffectType.Eradicant: FungicideEffect(1, 12.4)
            }
        )
        self.lowRisk = Fungicide(
            0.173,
            **{
                FungicideEffectType.Protectant: FungicideEffect(1, 2.1),
                FungicideEffectType.Eradicant: FungicideEffect(1, 2.1)
            }
        )
        self.psi = 1E-10
        self.phi = 0.13
        self.S0 = 42.34
        self._IS0 = 0
        self._IR0 = 0
        self.R0 = 0
        self.O0 = 0
        self.m = 0.1
        self.r = [0.147, 0.032]
        self.k = [26106, 2.461E8]
        self.seasonEndTime = 245  # Takes the simulation into September - which is traditional harvest time
        self.budBreakTime = 119
        self.floweringTime = 163
        self.veraisonTime = self.floweringTime + 7 * 7  # 6 - 8 weeks after flowering
        self.shootToppingTime = 173 + 1  # Shoot-topping should register the day after
        self.toppingLoss = {
            "S": 0.2,
            "E": 0.2,
            "I": 0.2,
            "R": 0.2,
            "O": 0.2
        }
        self.yieldStartTime = 163
        self.yieldEndTime = 193  # 30 days after flowering
        # self.sprayTimes = [119, 130, 141, 152, 163, 174, 185]
        # self.sprayTimes = [119, 133, 147, 161, 175]
        self.sprayTimes = [161, 175]
        self.alternatingSprayTimes = [161, 175]
        self.criticalYieldLoss = 0.03
        self.maxTime = self.seasonEndTime - self.budBreakTime
        self.rFreq = self.psi
        self._densityDependentInfection = True
        self._fungicideDecay = True
        self._seasonality = True
        self._latentPeriod = True

        self._setFittedParams()

    commonInputAttr = [
        "model", "gamma", "mu", "highRisk", "lowRisk", "psi", "phi",
        "S0", "ER0", "ES0", "IR0", "IS0", "R0", "r", "k", "criticalYieldLoss",
        "beta", "_densityDependentInfection", "_fungicideDecay", "_seasonality",
        "_latentPeriod", "maxTime", "highRiskDose", "lowRiskDose", "strategy"
    ]

    grapeOnlyInputAttr = [
        "O0", "m", "seasonEndTime", "budBreakTime", "floweringTime", "shootToppingTime",
        "toppingLoss", "sprayTimes"
    ]

    wheatOnlyInputAttr = [
        "nu", "PR0", "PS0", "GS32", "GS39", "GS61", "GS87", "nu"
    ]

    def hasSameInputAs(self, other, exceptFor = []) -> bool:
        compareOn = copy.deepcopy(Parameters.commonInputAttr)

        if self.model == Model.WheatSeptoria:
            compareOn += Parameters.wheatOnlyInputAttr
        else:
            compareOn += Parameters.grapeOnlyInputAttr

        # Remove anything we don't care about
        for exception in exceptFor:
            compareOn.remove(exception)

        return reduce(
            lambda left, el: left and getattr(self, el) == getattr(other, el),
            compareOn,
            True
        )

    def _clearAttributes(self):
        self.beta = None
        # TODO:

Parameters.betaMap.addModel(Model.WheatSeptoria)
Parameters.betaMap.addModel(Model.GrapePowderyMildew)

Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = True,
    seasonality = False,
    fungicideDecay = True,
    latentPeriod = True,
    value = [0.5593]
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = True,
    seasonality = False,
    fungicideDecay = True,
    latentPeriod = True,
    value = 0.0119
)
Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = False,
    seasonality = True,
    fungicideDecay = False,
    latentPeriod = True,
    value = [0.854, 0.2253]
)
Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = False,
    seasonality = True,
    fungicideDecay = True,
    latentPeriod = True,
    value = [0.9277, 0.1928]
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = True,
    seasonality = False,
    fungicideDecay = False,
    latentPeriod = False,
    value = 0.00484
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = False,
    seasonality = True,
    fungicideDecay = True,
    latentPeriod = False,
    value = 0.00548
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = False,
    seasonality = False,
    fungicideDecay = False,
    latentPeriod = False,
    value = 0.00464
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = False,
    seasonality = False,
    fungicideDecay = False,
    latentPeriod = True,
    value = 0.00972
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = False,
    seasonality = False,
    fungicideDecay = True,
    latentPeriod = True,
    value = 0.0109
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = True,
    seasonality = False,
    fungicideDecay = True,
    latentPeriod = False,
    value = 0.00535
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = True,
    seasonality = True,
    fungicideDecay = False,
    latentPeriod = True,
    value = 0.0135
)
Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = True,
    seasonality = True,
    fungicideDecay = True,
    latentPeriod = True,
    value = [1.605, 1.688]
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = False,
    seasonality = True,
    fungicideDecay = True,
    latentPeriod = True,
    value = 0.0127
)
Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = True,
    seasonality = False,
    fungicideDecay = False,
    latentPeriod = False,
    value = [0.2182]
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = True,
    seasonality = False,
    fungicideDecay = False,
    latentPeriod = True,
    value = 0.0103
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = True,
    seasonality = True,
    fungicideDecay = True,
    latentPeriod = False,
    value = 0.00688
)
Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = True,
    seasonality = True,
    fungicideDecay = True,
    latentPeriod = False,
    value = [0.4782, 0.913]
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = False,
    seasonality = True,
    fungicideDecay = False,
    latentPeriod = True,
    value = 0.0115
)
Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = False,
    seasonality = False,
    fungicideDecay = False,
    latentPeriod = True,
    value = [0.4554]
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = True,
    seasonality = True,
    fungicideDecay = False,
    latentPeriod = False,
    value = 0.00625
)
Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = False,
    seasonality = False,
    fungicideDecay = True,
    latentPeriod = False,
    value = [0.2015]
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = True,
    seasonality = True,
    fungicideDecay = True,
    latentPeriod = True,
    value = 1.56E-2
)
Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = True,
    seasonality = False,
    fungicideDecay = True,
    latentPeriod = False,
    value = [0.2336]
)
Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = False,
    seasonality = False,
    fungicideDecay = True,
    latentPeriod = True,
    value = [0.4569]
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = False,
    seasonality = True,
    fungicideDecay = False,
    latentPeriod = False,
    value = 0.0052
)
Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = True,
    seasonality = True,
    fungicideDecay = False,
    latentPeriod = True,
    value = [1.361, 1.637]
)
Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = False,
    seasonality = False,
    fungicideDecay = False,
    latentPeriod = False,
    value = [0.2025]
)
Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = False,
    seasonality = True,
    fungicideDecay = False,
    latentPeriod = False,
    value = [0.2711, 0.1487]
)
Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = False,
    seasonality = True,
    fungicideDecay = True,
    latentPeriod = False,
    value = [0.2817, 0.1376]
)
Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = True,
    seasonality = False,
    fungicideDecay = False,
    latentPeriod = True,
    value = [0.5046]
)
Parameters.betaMap.set(
    model = Model.WheatSeptoria,
    densityDependentInfection = False,
    seasonality = False,
    fungicideDecay = True,
    latentPeriod = False,
    value = 0.00497
)
Parameters.betaMap.set(
    model = Model.GrapePowderyMildew,
    densityDependentInfection = True,
    seasonality = True,
    fungicideDecay = False,
    latentPeriod = False,
    value = [0.4155, 0.893]
)
