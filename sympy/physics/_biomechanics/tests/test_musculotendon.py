"""Tests for objects from `sympy/physics/_biomechanics/musculotendon.py`."""


import pytest

import sympy as sm
import sympy.physics.mechanics as me
import sympy.physics._biomechanics as bme


@pytest.mark.parametrize(
    'musculotendon_class',
    [
        bme.Brockie2021Musculotendon,
        bme.DeGroote2016Musculotendon,
        bme.Millard2013Musculotendon,
    ]
)
class TestMusculotendonArguments:
    """Test correct handling of arguments for `MusculotendonBase`."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Instantiate name, and origin and insertion point fixtures."""
        self.name = 'muscle'
        self.origin = me.Point('origin')
        self.insertion = me.Point('insertion')

    @pytest.mark.parametrize('name', ['muscle', 'musculotendon'])
    def test_nonzero_length_strings_are_valid_names(
        self,
        musculotendon_class,
        name,
    ):
        """A `str` with nonzero length is a valid musculotendon name."""
        muscle = musculotendon_class(
            name,
            origin=self.origin,
            insertion=self.insertion,
        )
        assert muscle.name == name
        assert isinstance(muscle.name, str)

    @pytest.mark.parametrize('name', [0, sm.Symbol('name')])
    def test_non_string_name_raises_type_error(
        self,
        musculotendon_class,
        name,
    ):
        """A `TypeError` is raised if `name` is not a `str`."""
        with pytest.raises(TypeError):
            _ = musculotendon_class(
                name,
                origin=self.origin,
                insertion=self.insertion,
            )

    @pytest.mark.parametrize('name', [''])
    def test_invalid_string_name_raises_value_error(
        self,
        musculotendon_class,
        name,
    ):
        """A `ValueError` is raised if `name` is an invalid `str`."""
        with pytest.raises(ValueError):
            _ = musculotendon_class(
                name,
                origin=self.origin,
                insertion=self.insertion,
            )

    def test_point_objects_are_valid_origin_and_insertion(
        self,
        musculotendon_class,
    ):
        """`Point` objects passed to `origin` and `insertion` are valid."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
        )
        assert isinstance(muscle.origin, me.Point)
        assert muscle.origin == self.origin
        assert isinstance(muscle.insertion, me.Point)
        assert muscle.insertion == self.insertion

    @pytest.mark.parametrize('origin', ['origin', me.ReferenceFrame('N')])
    def test_origin_not_point_object_raises_type_error(
        self,
        musculotendon_class,
        origin,
    ):
        """A `TypeError` is raised if `origin` is not a `Point`."""
        with pytest.raises(TypeError):
            _ = musculotendon_class(
                self.name,
                origin=origin,
                insertion=self.insertion,
            )

    @pytest.mark.parametrize('insertion', ['insertion', me.ReferenceFrame('N')])
    def test_insertion_not_point_object_raises_type_error(
        self,
        musculotendon_class,
        insertion,
    ):
        """A `TypeError` is raised if `insertion` is not a `Point`."""
        with pytest.raises(TypeError):
            _ = musculotendon_class(
                self.name,
                origin=self.origin,
                insertion=insertion,
            )

    def test_optional_keyword_argument_defaults(self, musculotendon_class):
        """All parametrized symbol attributes default to `None` or value."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
        )
        assert muscle.optimal_fiber_length is None
        assert muscle.maximal_fiber_velocity == 10.0
        assert muscle.peak_isometric_force is None
        assert muscle.tendon_slack_length is None
        assert muscle.optimal_pennation_angle == 0.0
        assert muscle.fiber_damping_coefficient == 0.1

    def test_optimal_fiber_length_is_optional(self, musculotendon_class):
        """`None` is valid value for `optimal_fiber_length`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            optimal_fiber_length=None,
        )
        assert muscle.optimal_fiber_length is None

    @pytest.mark.parametrize(
        'optimal_fiber_length_value, optimal_fiber_length_expected',
        [
            (0.25, 0.25),
            (1, 1.0),
            ('0.25', 0.25),
            (sm.Float(0.25), 0.25),
        ]
    )
    def test_optimal_fiber_length_attribute_returns_float(
        self,
        musculotendon_class,
        optimal_fiber_length_value,
        optimal_fiber_length_expected,
    ):
        """`optimal_fiber_length` arguments are cast to `float`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            optimal_fiber_length=optimal_fiber_length_value,
        )
        assert isinstance(muscle.optimal_fiber_length, float)
        assert muscle.optimal_fiber_length == optimal_fiber_length_expected

    @pytest.mark.parametrize('optimal_fiber_length', [0.0, -0.25])
    def test_nonpositive_optimal_fiber_length_raises_value_error(
        self,
        musculotendon_class,
        optimal_fiber_length,
    ):
        """`ValueError` is raised for nonpositive `optimal_fiber_length`."""
        with pytest.raises(ValueError):
            _ = musculotendon_class(
                self.name,
                origin=self.origin,
                insertion=self.insertion,
                optimal_fiber_length=optimal_fiber_length,
            )

    def test_none_maximal_fiber_velocity_raises_type_error(
        self,
        musculotendon_class,
    ):
        """Passing `None` to `maximal_fiber_velocity` raises `TypeError`."""
        with pytest.raises(TypeError):
            _ = musculotendon_class(
                self.name,
                origin=self.origin,
                insertion=self.insertion,
                maximal_fiber_velocity=None,
            )

    @pytest.mark.parametrize(
        'maximal_fiber_velocity_value, maximal_fiber_velocity_expected',
        [
            (10.0, 10.0),
            (10, 10.0),
            ('10.0', 10.0),
            (sm.Float(10.0), 10.0),
            (sm.Integer(10), 10.0),
        ]
    )
    def test_maximal_fiber_velocity_attribute_returns_float(
        self,
        musculotendon_class,
        maximal_fiber_velocity_value,
        maximal_fiber_velocity_expected,
    ):
        """`maximal_fiber_velocity` arguments are cast to `float`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            maximal_fiber_velocity=maximal_fiber_velocity_value,
        )
        assert isinstance(muscle.maximal_fiber_velocity, float)
        assert muscle.maximal_fiber_velocity == maximal_fiber_velocity_expected

    @pytest.mark.parametrize('maximal_fiber_velocity', [0.0, -10.0])
    def test_nonpositive_maximal_fiber_velocity_raises_value_error(
        self,
        musculotendon_class,
        maximal_fiber_velocity,
    ):
        """`ValueError` is raised for nonpositive `maximal_fiber_velocity`."""
        with pytest.raises(ValueError):
            _ = musculotendon_class(
                self.name,
                origin=self.origin,
                insertion=self.insertion,
                maximal_fiber_velocity=maximal_fiber_velocity,
            )

    def test_peak_isometric_force_is_optional(self, musculotendon_class):
        """`None` is valid value for `peak_isometric_force`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            peak_isometric_force=None,
        )
        assert muscle.peak_isometric_force is None

    @pytest.mark.parametrize(
        'peak_isometric_force_value, peak_isometric_force_expected',
        [
            (1000.0, 1000.0),
            (1000, 1000.0),
            ('1000.0', 1000.0),
            (sm.Float(1000.0), 1000.0),
            (sm.Integer(1000), 1000.0)
        ]
    )
    def test_peak_isometric_force_attribute_returns_float(
        self,
        musculotendon_class,
        peak_isometric_force_value,
        peak_isometric_force_expected,
    ):
        """`peak_isometric_force` arguments are cast to `float`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            peak_isometric_force=peak_isometric_force_value,
        )
        assert isinstance(muscle.peak_isometric_force, float)
        assert muscle.peak_isometric_force == peak_isometric_force_expected

    @pytest.mark.parametrize('peak_isometric_force', [0.0, -1000.0])
    def test_nonpositive_peak_isometric_force_raises_value_error(
        self,
        musculotendon_class,
        peak_isometric_force,
    ):
        """`ValueError` is raised for nonpositive `peak_isometric_force`."""
        with pytest.raises(ValueError):
            _ = musculotendon_class(
                self.name,
                origin=self.origin,
                insertion=self.insertion,
                peak_isometric_force=peak_isometric_force,
            )

    def test_tendon_slack_length_is_optional(self, musculotendon_class):
        """`None` is valid value for `tendon_slack_length`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            tendon_slack_length=None,
        )
        assert muscle.tendon_slack_length is None

    @pytest.mark.parametrize(
        'tendon_slack_length_value, tendon_slack_length_expected',
        [
            (0.05, 0.05),
            (1, 1.0),
            ('0.05', 0.05),
            (sm.Float(0.05), 0.05),
        ]
    )
    def test_tendon_slack_length_attribute_returns_float(
        self,
        musculotendon_class,
        tendon_slack_length_value,
        tendon_slack_length_expected,
    ):
        """`tendon_slack_length` arguments are cast to `float`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            tendon_slack_length=tendon_slack_length_value,
        )
        assert isinstance(muscle.tendon_slack_length, float)
        assert muscle.tendon_slack_length == tendon_slack_length_expected

    @pytest.mark.parametrize('tendon_slack_length', [0.0, -0.05])
    def test_nonpositive_tendon_slack_length_raises_value_error(
        self,
        musculotendon_class,
        tendon_slack_length,
    ):
        """`ValueError` is raised for nonpositive `tendon_slack_length`."""
        with pytest.raises(ValueError):
            _ = musculotendon_class(
                self.name,
                origin=self.origin,
                insertion=self.insertion,
                tendon_slack_length=tendon_slack_length,
            )

    def test_none_optimal_pennation_angle_raises_type_error(
        self,
        musculotendon_class,
    ):
        """Passing `None` to `optimal_pennation_angle` raises `TypeError`."""
        with pytest.raises(TypeError):
            _ = musculotendon_class(
                self.name,
                origin=self.origin,
                insertion=self.insertion,
                optimal_pennation_angle=None,
            )

    @pytest.mark.parametrize(
        'optimal_pennation_angle_value, optimal_pennation_angle_expected',
        [
            (0.25, 0.25),
            (1, 1.0),
            ('0.25', 0.25),
            (sm.Float(0.25), 0.25),
        ]
    )
    def test_optimal_pennation_angle_attribute_returns_float(
        self,
        musculotendon_class,
        optimal_pennation_angle_value,
        optimal_pennation_angle_expected,
    ):
        """`optimal_pennation_angle` arguments are cast to `float`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            optimal_pennation_angle=optimal_pennation_angle_value,
        )
        assert isinstance(muscle.optimal_pennation_angle, float)
        assert muscle.optimal_pennation_angle == optimal_pennation_angle_expected

    @pytest.mark.parametrize('optimal_pennation_angle', [-0.25])
    def test_nonpositive_optimal_pennation_angle_raises_value_error(
        self,
        musculotendon_class,
        optimal_pennation_angle,
    ):
        """`ValueError` is raised for nonpositive `optimal_pennation_angle`."""
        with pytest.raises(ValueError):
            _ = musculotendon_class(
                self.name,
                origin=self.origin,
                insertion=self.insertion,
                optimal_pennation_angle=optimal_pennation_angle,
            )

    def test_none_fiber_damping_coefficient_raises_type_error(
        self,
        musculotendon_class,
    ):
        """Passing `None` to `fiber_damping_coefficient` raises `TypeError`."""
        with pytest.raises(TypeError):
            _ = musculotendon_class(
                self.name,
                origin=self.origin,
                insertion=self.insertion,
                fiber_damping_coefficient=None,
            )

    @pytest.mark.parametrize(
        'fiber_damping_coefficient_value, fiber_damping_coefficient_expected',
        [
            (0.1, 0.1),
            (1, 1.0),
            ('0.1', 0.1),
            (sm.Float(0.1), 0.1),
        ]
    )
    def test_fiber_damping_coefficient_attribute_returns_float(
        self,
        musculotendon_class,
        fiber_damping_coefficient_value,
        fiber_damping_coefficient_expected,
    ):
        """`fiber_damping_coefficient` arguments are numbers and cast to `float`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            fiber_damping_coefficient=fiber_damping_coefficient_value,
        )
        assert isinstance(muscle.fiber_damping_coefficient, float)
        assert muscle.fiber_damping_coefficient == fiber_damping_coefficient_expected

    @pytest.mark.parametrize('fiber_damping_coefficient', [-0.1])
    def test_nonpositive_fiber_damping_coefficient_raises_value_error(
        self,
        musculotendon_class,
        fiber_damping_coefficient,
    ):
        """`ValueError` is raised for nonpositive `fiber_damping_coefficient`."""
        with pytest.raises(ValueError):
            _ = musculotendon_class(
                self.name,
                origin=self.origin,
                insertion=self.insertion,
                fiber_damping_coefficient=fiber_damping_coefficient,
            )


@pytest.mark.parametrize(
    'musculotendon_class',
    [
        bme.Brockie2021Musculotendon,
        bme.DeGroote2016Musculotendon,
        bme.Millard2013Musculotendon,
    ]
)
class TestMusculotendonAttributes:
    """Tests for the attributes of `MusculotendonBase`."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Instantiate name, and origin and insertion point fixtures."""
        self.name = 'muscle'
        self.origin = me.Point('origin')
        self.insertion = me.Point('insertion')

    @pytest.mark.parametrize(
        'attribute_name, attribute_symbol',
        [
            ('l_M_opt', sm.Symbol('l_M_opt_muscle')),
            ('v_M_max', sm.Symbol('v_M_max_muscle')),
            ('F_M_max', sm.Symbol('F_M_max_muscle')),
            ('l_T_slack', sm.Symbol('l_T_slack_muscle')),
            ('alpha_opt', sm.Symbol('alpha_opt_muscle')),
            ('beta', sm.Symbol('beta_muscle')),
        ]
    )
    def test_has_symbol_as_attribute(
        self,
        musculotendon_class,
        attribute_name,
        attribute_symbol,
    ):
        """Attributes with `Symbol` instances exist for all expected names.

        Notes
        =====
        The symbols created are suffixed with the `name` associated with the
        musculotendon, which is `muscle` in this instance, separated from the
        symbol identifier by an underscore.

        """
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
        )
        assert hasattr(muscle, attribute_name)
        assert getattr(muscle, attribute_name) == attribute_symbol

    @pytest.mark.parametrize(
        'attribute_name, attribute_symbol',
        [
            ('l_MT', me.dynamicsymbols('l_MT_muscle')),
            ('v_MT', me.dynamicsymbols('v_MT_muscle')),
            ('l_T', me.dynamicsymbols('l_T_muscle')),
            ('v_T', me.dynamicsymbols('v_T_muscle')),
            ('l_M', me.dynamicsymbols('l_M_muscle')),
            ('v_M', me.dynamicsymbols('v_M_muscle')),
            ('l_T_tilde', me.dynamicsymbols('l_T_tilde_muscle')),
            ('v_T_tilde', me.dynamicsymbols('v_T_tilde_muscle')),
            ('l_M_tilde', me.dynamicsymbols('l_M_tilde_muscle')),
            ('v_M_tilde', me.dynamicsymbols('v_M_tilde_muscle')),
            ('F_T', me.dynamicsymbols('F_T_muscle')),
            ('F_M', me.dynamicsymbols('F_M_muscle')),
            ('F_T_tilde', me.dynamicsymbols('F_T_tilde_muscle')),
            ('F_M_tilde', me.dynamicsymbols('F_M_tilde_muscle')),
            ('fl_T', me.dynamicsymbols('fl_T_muscle')),
            ('fl_M_pas', me.dynamicsymbols('fl_M_pas_muscle')),
            ('fl_M_act', me.dynamicsymbols('fl_M_act_muscle')),
            ('fv_M', me.dynamicsymbols('fv_M_muscle')),
            ('cos_alpha', me.dynamicsymbols('cos_alpha_muscle')),
            ('dF_T_tilde_dt', me.dynamicsymbols('dF_T_tilde_dt_muscle')),
            ('dl_M_tilde_dt', me.dynamicsymbols('dl_M_tilde_dt_muscle')),
            ('F_orig_x', me.dynamicsymbols('F_orig_x_muscle')),
            ('F_orig_y', me.dynamicsymbols('F_orig_y_muscle')),
            ('F_orig_z', me.dynamicsymbols('F_orig_z_muscle')),
            ('F_insr_x', me.dynamicsymbols('F_insr_x_muscle')),
            ('F_insr_y', me.dynamicsymbols('F_insr_y_muscle')),
            ('F_insr_z', me.dynamicsymbols('F_insr_z_muscle')),
        ]
    )
    def test_has_dynamic_symbol_as_attribute(
        self,
        musculotendon_class,
        attribute_name,
        attribute_symbol,
    ):
        """Attributes with `Symbol` instances exist for all expected names.

        Notes
        =====
        The symbols created are suffixed with the `name` associated with the
        musculotendon, which is `muscle` in this instance, separated from the
        symbol identifier by an underscore.

        """
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
        )
        assert hasattr(muscle, attribute_name)
        assert getattr(muscle, attribute_name) == attribute_symbol

    def test_has_symbol_to_constant_mapping(self, musculotendon_class):
        """Has attribute that maps symbol attributes to constant attributes."""
        OPTIMAL_FIBER_LENGTH = 0.25
        MAXIMAL_FIBER_VELOCITY = 9.0
        PEAK_ISOMETRIC_FORCE = 1000.0
        TENDON_SLACK_LENGTH = 0.05
        OPTIMAL_PENNATION_ANGLE = 0.05
        FIBER_DAMPING_COEFFICIENT = 0.2
        muscle = musculotendon_class(
            'muscle',
            origin=self.origin,
            insertion=self.insertion,
            optimal_fiber_length=OPTIMAL_FIBER_LENGTH,
            maximal_fiber_velocity=MAXIMAL_FIBER_VELOCITY,
            peak_isometric_force=PEAK_ISOMETRIC_FORCE,
            tendon_slack_length=TENDON_SLACK_LENGTH,
            optimal_pennation_angle=OPTIMAL_PENNATION_ANGLE,
            fiber_damping_coefficient=FIBER_DAMPING_COEFFICIENT,
        )
        expected = {
            sm.Symbol('l_M_opt_muscle') :OPTIMAL_FIBER_LENGTH,
            sm.Symbol('v_M_max_muscle') :MAXIMAL_FIBER_VELOCITY,
            sm.Symbol('F_M_max_muscle') :PEAK_ISOMETRIC_FORCE,
            sm.Symbol('l_T_slack_muscle') :TENDON_SLACK_LENGTH,
            sm.Symbol('alpha_opt_muscle') :OPTIMAL_PENNATION_ANGLE,
            sm.Symbol('beta_muscle') :FIBER_DAMPING_COEFFICIENT,
        }
        assert hasattr(muscle, 'symbol_to_constant_mapping')
        assert muscle.symbol_to_constant_mapping() == expected


class TestMusculotendonFactoryFunction:
    """Tests for the `Musculotendon` factory function."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Instantiate origin and insertion point fixtures."""
        self.origin = me.Point('origin')
        self.insertion = me.Point('insertion')

    @pytest.mark.parametrize(
        'model_identifier, musculotendon_class',
        [
            ('Brockie', bme.Brockie2021Musculotendon),
            ('Brockie2021', bme.Brockie2021Musculotendon),
            ('DeGroote', bme.DeGroote2016Musculotendon),
            ('DeGroote2016', bme.DeGroote2016Musculotendon),
            ('Millard', bme.Millard2013Musculotendon),
            ('Millard2013', bme.Millard2013Musculotendon),
        ]
    )
    def test_valid_identifier_returns_musclotendon_instance(
        self,
        model_identifier,
        musculotendon_class,
    ):
        """A correct identifier will successfully return an instance."""
        muscle = bme.Musculotendon(
            'muscle',
            model_identifier,
            origin=self.origin,
            insertion=self.insertion,
        )
        assert isinstance(muscle, musculotendon_class)

    @pytest.mark.parametrize('model_identifier', ['invalid_identifier'])
    def test_invalid_identifier_raises_value_error(
        self,
        model_identifier,
    ):
        """An invalid identifier will cause a `ValueError` to be raised."""
        with pytest.raises(ValueError):
            _ = bme.Musculotendon(
                'muscle',
                model_identifier,
                origin=self.origin,
                insertion=self.insertion,
            )
