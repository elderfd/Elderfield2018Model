#!usr/bin/python


class MetaConst(type):
    def __getattr__(cls, key):
        return cls[key]

    def __setattr__(cls, key, value):
        raise TypeError


class Const(object, metaclass = MetaConst):
    def __getattr__(self, name):
        return self[name]

    def __setattr__(self, name, value):
        raise TypeError


class SprayStrategies(Const):
    Mixture = "mixture"
    AltHiLo = "alternationHiLo"
    AltLoHi = "alternationLoHi"
    NoSpray = "noSpray"


class Model(Const):
    WheatSeptoria = "wheatSeptoria"
    GrapePowderyMildew = "grapePowderyMildew"


class StateIndices(Const):
    S, ER, ES, IR, IS, R, PR, PS, O, high, low = range(11)


class FungicideEffectType(Const):
    Eradicant = "Eradicant"
    Protectant = "Protectant"
