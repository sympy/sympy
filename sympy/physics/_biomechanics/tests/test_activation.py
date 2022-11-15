"""Tests for objects from `sympy/physics/_biomechanics/activation.py`."""


import pytest

import sympy as sm
import sympy.physics.mechanics as me
import sympy.physics._biomechanics as bme


EXPECTED_ATTRIBUTE_MISSING_SENTINEL = object()


@pytest.mark.parametrize(
    'activation_dynamics_class',
    [
        bme.DeGroote2016ActivationDynamics,
        bme.ZerothOrderActivationDynamics,
    ]
)
class TestActivationDynamicsArguments:
    """Test correct handling of arguments for `ActivationDynamicsBase`."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Instantiate name."""
        self.name = 'activation'

    @pytest.mark.parametrize('name', ['act', 'activation'])
    def test_nonzero_length_strings_are_valid_names(
        self,
        activation_dynamics_class,
        name,
    ):
        """A `str` with nonzero length is a valid musculotendon name."""
        activation = activation_dynamics_class(name)
        assert activation.name == name
        assert isinstance(activation.name, str)

    @pytest.mark.parametrize('name', [0, sm.Symbol('name')])
    def test_non_string_name_raises_type_error(
        self,
        activation_dynamics_class,
        name,
    ):
        """A `TypeError` is raised if `name` is not a `str`."""
        with pytest.raises(TypeError):
            _ = activation_dynamics_class(name)

    @pytest.mark.parametrize('name', [''])
    def test_invalid_string_name_raises_value_error(
        self,
        activation_dynamics_class,
        name,
    ):
        """A `ValueError` is raised if `name` is an invalid `str`."""
        with pytest.raises(ValueError):
            _ = activation_dynamics_class(name)


@pytest.mark.parametrize(
    'activation_dynamics_class, attributes',
    [
        (
            bme.DeGroote2016ActivationDynamics,
            {
                "order": 1,
            },
        ),
        (
            bme.ZerothOrderActivationDynamics,
            {
                "order": 0,
            },
        ),
    ]
)
class TestActivationDynamicsAttributes:
    """Test attributes and properties of `ActivationDynamicsBase`."""

    @pytest.fixture(autouse=True)
    def setup(self, activation_dynamics_class, attributes):
        """Create basic instance."""
        self.instance = activation_dynamics_class('activation')
        self.attributes = attributes

    def test_order_class_attribute(self):
        """Instances have an `int` `order` attribute."""
        expected = self.attributes.get(
            'order',
            EXPECTED_ATTRIBUTE_MISSING_SENTINEL,
        )
        assert self.instance.order == expected
        assert isinstance(self.instance.order, int)

    @pytest.mark.parametrize(
        'attribute_name, attribute_symbol',
        [
            ('e', me.dynamicsymbols('e_activation')),
            ('a', me.dynamicsymbols('a_activation')),
        ]
    )
    def test_has_symbol_as_attribute(
        self,
        attribute_name,
        attribute_symbol,
    ):
        """Attributes with `Symbol` instances exist for all expected names.

        These symbol instances are actually `dynamicsymbols` as they need to be
        functions of time.

        Notes
        =====
        The dynamic symbols created are suffixed with the `name` associated with
        the activation dynamics, which is `activation` in this instance,
        separated from the symbol identifier by an underscore.

        """
        assert hasattr(self.instance, attribute_name)
        assert getattr(self.instance, attribute_name) == attribute_symbol


class TestZerothOrderActivationDynamics:
    """Tests specific to the `ZerothOrderActivationDynamics` class."""

    def test_has_symbol_to_constant_mapping(self):
        """Has attribute that maps symbol attributes to constant attributes."""
        activation = bme.ZerothOrderActivationDynamics('activation')
        expected = {}
        assert hasattr(activation, 'symbol_to_constant_mapping')
        assert activation.symbol_to_constant_mapping() == expected


class TestDeGroote2016ActivationDynamics:
    """Tests specific to the `DeGroote2016ActivationDynamics` class."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Create basic instance."""
        self.name = 'activation'
        self.instance = bme.DeGroote2016ActivationDynamics(self.name)

    def test_optional_keyword_argument_default_values(self):
        """All optional keyword arguments default to expected values."""
        instance = bme.DeGroote2016ActivationDynamics('activation')
        assert instance.activation_time_constant == 0.015
        assert instance.deactivation_time_constant == 0.060

    @pytest.mark.parametrize(
        'activation_time_constant_value, activation_time_constant_expected',
        [
            (0.020, 0.020),
            (1, 1.0),
            ('0.020', 0.020),
            (sm.Float(0.020), 0.020),
        ]
    )
    def test_activation_time_constant_attribute_returns_float(
        self,
        activation_time_constant_value,
        activation_time_constant_expected,
    ):
        """`activation_time_constant` arguments are cast to `float`."""
        instance = bme.DeGroote2016ActivationDynamics(
            self.name,
            activation_time_constant=activation_time_constant_value,
        )
        assert isinstance(instance.activation_time_constant, float)
        assert instance.activation_time_constant == activation_time_constant_expected

    @pytest.mark.parametrize('activation_time_constant', [0.0, -0.25])
    def test_nonpositive_activation_time_constant_raises_value_error(
        self,
        activation_time_constant,
    ):
        """`ValueError` is raised for nonpositive `activation_time_constant`."""
        with pytest.raises(ValueError):
            _ = bme.DeGroote2016ActivationDynamics(
                self.name,
                activation_time_constant=activation_time_constant,
            )

    @pytest.mark.parametrize(
        'attribute_name, attribute_symbol',
        [
            ('tau_a', sm.Symbol('tau_a_activation')),
            ('tau_d', sm.Symbol('tau_d_activation')),
        ]
    )
    def test_has_symbol_as_attribute(
        self,
        attribute_name,
        attribute_symbol,
    ):
        """Attributes with `Symbol` instances exist for all expected names.

        Notes
        =====
        The dynamic symbols created are suffixed with the `name` associated with
        the activation dynamics, which is `activation` in this instance,
        separated from the symbol identifier by an underscore.

        """
        assert hasattr(self.instance, attribute_name)
        assert getattr(self.instance, attribute_name) == attribute_symbol

    def test_has_symbol_to_constant_mapping(self):
        """Has attribute that maps symbol attributes to constant attributes."""
        ACTIVATION_TIME_CONSTANT = 0.025
        DEACTIVATION_TIME_CONSTANT = 0.070
        activation = bme.DeGroote2016ActivationDynamics(
            'activation',
            activation_time_constant=ACTIVATION_TIME_CONSTANT,
            deactivation_time_constant=DEACTIVATION_TIME_CONSTANT,
        )
        expected = {
            sm.Symbol('tau_a_activation'): ACTIVATION_TIME_CONSTANT,
            sm.Symbol('tau_d_activation'): DEACTIVATION_TIME_CONSTANT,
        }
        assert hasattr(activation, 'symbol_to_constant_mapping')
        assert activation.symbol_to_constant_mapping() == expected


class TestActivationDynamicsFactoryFunction:
    """Tests for the `ActivationDynamics` factory function."""

    @pytest.mark.parametrize(
        'model_identifier, activation_dynamics_class',
        [
            (0, bme.ZerothOrderActivationDynamics),
            ('Zero', bme.ZerothOrderActivationDynamics),
            ('Zeroth', bme.ZerothOrderActivationDynamics),
            (1, bme.DeGroote2016ActivationDynamics),
            ('DeGroote', bme.DeGroote2016ActivationDynamics),
            ('DeGroote2016', bme.DeGroote2016ActivationDynamics),
        ]
    )
    def test_valid_identifier_returns_activation_dynamics_instance(
        self,
        model_identifier,
        activation_dynamics_class,
    ):
        """A correct identifier will successfully return an instance."""
        activation = bme.ActivationDynamics(
            'activation',
            model_identifier,
        )
        assert isinstance(activation, activation_dynamics_class)

    @pytest.mark.parametrize('model_identifier', [2, 'invalid_identifier'])
    def test_invalid_identifier_raises_value_error(
        self,
        model_identifier,
    ):
        """An invalid identifier will cause a `ValueError` to be raised."""
        with pytest.raises(ValueError):
            _ = bme.ActivationDynamics(
                'activation',
                model_identifier,
            )
