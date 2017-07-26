from .constants import Model
from .stepChange import StepChange, StepChangeKind
from .constants import StateIndices
from operator import attrgetter
from . import types
import typing
from .parameters import Parameters


class FungicideTreatment(object):
    def __init__(
        self,
        highRiskTimes: typing.List[types.Number],
        highRiskDoses: typing.List[types.Number],
        lowRiskTimes: typing.List[types.Number],
        lowRiskDoses: typing.List[types.Number]
    ) -> None:
        self.highRiskTimes = highRiskTimes
        self.highRiskDoses = highRiskDoses
        self.lowRiskTimes = lowRiskTimes
        self.lowRiskDoses = lowRiskDoses

    def createStepChanges(self, params: Parameters) -> typing.List[StepChange]:
        if params.model != Model.WheatSeptoria:
            raise NotImplementedError("Object-based treatments only implemented for wheat model.")

        t = 0
        changes = []

        if not params.fungicideDecay:
            highRiskLifetime = 1 / params.highRisk.decayRate
            lowRiskLifetime = 1 / params.lowRisk.decayRate

            maxLifeTime = params.GS39 - params.GS32
            largestLifeTime = max(highRiskLifetime, lowRiskLifetime)

            if largestLifeTime >= maxLifeTime:
                multiplier = maxLifeTime / largestLifeTime
                highRiskLifetime *= multiplier
                lowRiskLifetime *= multiplier
                highRiskLifetime = round(highRiskLifetime, 6)
                lowRiskLifetime = round(lowRiskLifetime, 6)

        while t < params.maxTime:
            # Do start of season initialisation
            changes.append(StepChange(
                StateIndices.S, t, params.S0, StepChangeKind.SET, params
            ))
            changes.append(StepChange(
                StateIndices.ER, t, attrgetter("ER0"), StepChangeKind.SET, params
            ))
            changes.append(StepChange(
                StateIndices.ES, t, attrgetter("ES0"), StepChangeKind.SET, params
            ))
            changes.append(StepChange(
                StateIndices.IR, t, attrgetter("IR0"), StepChangeKind.SET, params
            ))
            changes.append(StepChange(
                StateIndices.IS, t, attrgetter("IS0"), StepChangeKind.SET, params
            ))
            changes.append(StepChange(
                StateIndices.R, t, params.R0, StepChangeKind.SET, params
            ))
            changes.append(StepChange(
                StateIndices.PR, t, attrgetter("PR0"), StepChangeKind.SET, params
            ))
            changes.append(StepChange(
                StateIndices.PS, t, attrgetter("PS0"), StepChangeKind.SET, params
            ))

            for sprayTime, dose in zip(self.lowRiskTimes, self.lowRiskDoses):
                changes.append(StepChange(
                    StateIndices.low, t + sprayTime, dose, StepChangeKind.ADD, params
                ))
                if not params.fungicideDecay:
                    changes.append(StepChange(
                        StateIndices.low, t + sprayTime + lowRiskLifetime, 0, StepChangeKind.SET, params
                    ))
            for sprayTime, dose in zip(self.highRiskTimes, self.highRiskDoses):
                changes.append(StepChange(
                    StateIndices.high, t + sprayTime, dose, StepChangeKind.ADD, params
                ))
                if not params.fungicideDecay:
                    changes.append(StepChange(
                        StateIndices.high, t + sprayTime + highRiskLifetime, 0, StepChangeKind.SET, params
                    ))

            # To next season!
            t += params.GS87 - params.seasonStartTime

        return changes
