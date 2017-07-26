#!usr/bin/python
from math import exp
from . import types
import typing


class FungicideEffect(object):
    def __init__(self, maxEffect: types.Number, curvatureFactor: types.Number) -> None:
        self.maxEffect = maxEffect
        self.curvatureFactor = curvatureFactor

    def effect(self, dose: types.Number) -> types.Number:
        return self.maxEffect * (1 - exp(-self.curvatureFactor * dose))


class Fungicide(object):
    def __init__(self, decayRate: float, **effects: FungicideEffect) -> None:
        self.decayRate = decayRate
        self.effects = {k: effects[k] for k in effects}

    def effect(self, dose: types.Number, effect: str) -> types.Number:
        if effect not in self.effects:
            return 0
        else:
            return self.effects[effect].effect(dose)

    @staticmethod
    def combineEffects(fungicides: typing.List["Fungicide"], doses: typing.List[types.Number], effect: str):
        complement = 1

        for fungicide, dose in zip(fungicides, doses):
            complement *= (1 - fungicide.effect(dose, effect))

        return 1 - complement

    def __eq__(self, other: object) -> bool:
        return self.__dict__ == other.__dict__
