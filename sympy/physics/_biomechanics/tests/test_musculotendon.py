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
    """Test correct handling of arguments for `MusculotendonABC`."""

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

    def test_none_optimal_fiber_length_instantiates_default_symbol(
        self,
        musculotendon_class,
    ):
        """Passing `None` to `optimal_fiber_length` creates `Symbol`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            optimal_fiber_length=None,
        )
        assert isinstance(muscle.optimal_fiber_length, sm.Symbol)
        assert muscle.optimal_fiber_length == sm.Symbol('l_M_opt_muscle')

    @pytest.mark.parametrize('optimal_fiber_length', [(0.25), (sm.Float(0.25))])
    def test_numeric_optimal_fiber_length_instantiates_number_symbol(
        self,
        musculotendon_class,
        optimal_fiber_length,
    ):
        """Passing a numeric value to `optimal_fiber_length` creates `Number`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            optimal_fiber_length=optimal_fiber_length,
        )
        assert isinstance(muscle.optimal_fiber_length, sm.Number)
        assert muscle.optimal_fiber_length == sm.Number(optimal_fiber_length)

    def test_none_maximal_fiber_velocity_instantiates_default_symbol(
        self,
        musculotendon_class,
    ):
        """Passing `None` to `maximal_fiber_velocity` creates `Symbol`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            maximal_fiber_velocity=None,
        )
        assert isinstance(muscle.maximal_fiber_velocity, sm.Symbol)
        assert muscle.maximal_fiber_velocity == sm.Symbol('v_M_max_muscle')

    @pytest.mark.parametrize('maximal_fiber_velocity', [(10.0), (sm.Float(10.0))])
    def test_numeric_maximal_fiber_velocity_instantiates_number_symbol(
        self,
        musculotendon_class,
        maximal_fiber_velocity,
    ):
        """Passing a numeric value to `maximal_fiber_velocity` creates `Number`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            maximal_fiber_velocity=maximal_fiber_velocity,
        )
        assert isinstance(muscle.maximal_fiber_velocity, sm.Number)
        assert muscle.maximal_fiber_velocity == sm.Number(maximal_fiber_velocity)

    def test_none_peak_isometric_force_instantiates_default_symbol(
        self,
        musculotendon_class,
    ):
        """Passing `None` to `peak_isometric_force` creates `Symbol`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            peak_isometric_force=None,
        )
        assert isinstance(muscle.peak_isometric_force, sm.Symbol)
        assert muscle.peak_isometric_force == sm.Symbol('F_M_max_muscle')

    @pytest.mark.parametrize('peak_isometric_force', [(1000.0), (sm.Float(1000.0))])
    def test_numeric_peak_isometric_force_instantiates_number_symbol(
        self,
        musculotendon_class,
        peak_isometric_force,
    ):
        """Passing a numeric value to `peak_isometric_force` creates `Number`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            peak_isometric_force=peak_isometric_force,
        )
        assert isinstance(muscle.peak_isometric_force, sm.Number)
        assert muscle.peak_isometric_force == sm.Number(peak_isometric_force)

    def test_none_tendon_slack_length_instantiates_default_symbol(
        self,
        musculotendon_class,
    ):
        """Passing `None` to `tendon_slack_length` creates `Symbol`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            tendon_slack_length=None,
        )
        assert isinstance(muscle.tendon_slack_length, sm.Symbol)
        assert muscle.tendon_slack_length == sm.Symbol('l_T_slack_muscle')

    @pytest.mark.parametrize('tendon_slack_length', [(0.05), (sm.Float(0.05))])
    def test_numeric_tendon_slack_length_instantiates_number_symbol(
        self,
        musculotendon_class,
        tendon_slack_length,
    ):
        """Passing a numeric value to `tendon_slack_length` creates `Number`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            tendon_slack_length=tendon_slack_length,
        )
        assert isinstance(muscle.tendon_slack_length, sm.Number)
        assert muscle.tendon_slack_length == sm.Number(tendon_slack_length)

    def test_none_optimal_pennation_angle_instantiates_default_symbol(
        self,
        musculotendon_class,
    ):
        """Passing `None` to `optimal_pennation_angle` creates `Symbol`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            optimal_pennation_angle=None,
        )
        assert isinstance(muscle.optimal_pennation_angle, sm.Symbol)
        assert muscle.optimal_pennation_angle == sm.Symbol('alpha_opt_muscle')

    @pytest.mark.parametrize('optimal_pennation_angle', [(0.0), (sm.Float(0.0))])
    def test_numeric_optimal_pennation_angle_instantiates_number_symbol(
        self,
        musculotendon_class,
        optimal_pennation_angle,
    ):
        """Passing a numeric value to `optimal_pennation_angle` creates `Number`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            optimal_pennation_angle=optimal_pennation_angle,
        )
        assert isinstance(muscle.optimal_pennation_angle, sm.Number)
        assert muscle.optimal_pennation_angle == sm.Number(optimal_pennation_angle)

    def test_none_fiber_damping_coefficient_instantiates_default_symbol(
        self,
        musculotendon_class,
    ):
        """Passing `None` to `fiber_damping_coefficient` creates `Symbol`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            fiber_damping_coefficient=None,
        )
        assert isinstance(muscle.fiber_damping_coefficient, sm.Symbol)
        assert muscle.fiber_damping_coefficient == sm.Symbol('beta_muscle')

    @pytest.mark.parametrize('fiber_damping_coefficient', [(0.1), (sm.Float(0.1))])
    def test_numeric_fiber_damping_coefficient_instantiates_number_symbol(
        self,
        musculotendon_class,
        fiber_damping_coefficient,
    ):
        """Passing a numeric value to `fiber_damping_coefficient` creates `Number`."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
            fiber_damping_coefficient=fiber_damping_coefficient,
        )
        assert isinstance(muscle.fiber_damping_coefficient, sm.Number)
        assert muscle.fiber_damping_coefficient == sm.Number(fiber_damping_coefficient)


@pytest.mark.parametrize(
    'musculotendon_class',
    [
        bme.Brockie2021Musculotendon,
        bme.DeGroote2016Musculotendon,
        bme.Millard2013Musculotendon,
    ]
)
class TestMusculotendonSymbolicAttributes:
    """Tests for the symbolic attributes of `MusculotendonABC`."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Instantiate name, and origin and insertion point fixtures."""
        self.name = 'muscle'
        self.origin = me.Point('origin')
        self.insertion = me.Point('insertion')

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
    def test_has_symbol_as_attribute(
        self,
        musculotendon_class,
        attribute_name,
        attribute_symbol,
    ):
        """Attributes with `Symbol` instances exist for all expected names."""
        muscle = musculotendon_class(
            self.name,
            origin=self.origin,
            insertion=self.insertion,
        )
        assert hasattr(muscle, attribute_name)
        assert getattr(muscle, attribute_name) == attribute_symbol


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

    @pytest.mark.parametrize('model_identifier', ['invalid_name'])
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
