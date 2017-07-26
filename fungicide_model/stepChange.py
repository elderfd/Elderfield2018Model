#!usr/bin/python
from .constants import Const
from operator import attrgetter


class Change(object):
    def __init__(self, amount, kind, context):
        self.amount = amount
        self.kind = kind
        self.context = context

    def value(self):
        if isinstance(self.amount, attrgetter):
            return self.amount(self.context)
        elif callable(self.amount):
            return self.amount()
        else:
            return self.amount

    def apply(self, x):
        if self.kind == StepChangeKind.ADD:
            return self._add(x)
        elif self.kind == StepChangeKind.MULTIPLY:
            return self._multiply(x)
        else:
            return self._set()

    def _add(self, x):
        return x + self.value()

    def _set(self):
        return self.value()

    def _multiply(self, x):
        return self.value() * x


class StepChangeKind(Const):
    SET = 0
    ADD = 1
    MULTIPLY = 2


class StepChange(object):
    def __init__(self, index, time, amount, kind, context):
        self.index = index
        self.time = time
        self.change = Change(amount, kind, context)

    def apply(self, state):
        newState = state

        newState[self.index] = self.change.apply(state[self.index])

        return newState

    def __repr__(self):
        string = "Index: " + str(self.index) + ", Time: " + str(self.time)
        string += ", Change: "

        if self.change.apply(0) != self.change.apply(1):
            string += "+"

        string += str(self.change.apply(0))

        return string
