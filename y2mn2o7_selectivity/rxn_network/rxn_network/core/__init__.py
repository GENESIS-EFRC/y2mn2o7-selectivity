" Core interfaces for RXN Network"
from abc import ABCMeta, abstractmethod, abstractproperty
from monty.json import MSONable
from rxn_network.core.reaction import Reaction


class Calculator(MSONable, metaclass=ABCMeta):
    " Base definition for a property calculator "

    @abstractmethod
    def calculate(self, rxn: Reaction) -> float:
        " Evaluates the specified property of a reaction"


class CostFunction(MSONable, metaclass=ABCMeta):
    " Base definition for a cost function "

    @abstractmethod
    def evaluate(self, rxn: Reaction) -> float:
        " Evaluates the total cost function on a reaction "

