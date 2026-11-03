"""
Levenberg-Marquardt fitting engine.

Wrapper around scipy.optimize.leastsq with internal bound transformations.
"""

from typing import Callable, List, Dict, Optional, Tuple, Any
import numpy as np
from scipy.optimize import leastsq, OptimizeResult


class FitParameter:
    """
    Single fit parameter with bounds and constraints.

    Attributes
    ----------
    name : str
        Parameter name (e.g., 'a', 'b', 'c', 'd', 'e')
    value : float
        Current value
    min : float
        Lower bound
    max : float
        Upper bound
    fixed : bool
        If True, parameter is not varied during fit
    flags : tuple
        List of parameter flags
    """

    def __init__(
        self,
        name: str,
        value: float,
        min: float = -np.inf,
        max: float = np.inf,
        fixed: bool = False,
        flags: Tuple = (),
    ):
        self.name = name
        self.value = value
        self.min = min
        self.max = max
        self.fixed = fixed
        self.flags = flags

    def __repr__(self) -> str:
        status = "fixed" if self.fixed else "free"
        return f"FitParameter({self.name}={self.value:.4g}, [{self.min}, {self.max}], {status}, {self.flags})"


class Parameters:
    """
    Collection of fit parameters for CMWP.

    Manages all parameters, bounds, and conversions between
    dictionary and array representations for scipy.
    """

    def __init__(self):
        self.parameters: Dict[str, FitParameter] = {}

    def add(
        self,
        name: str,
        value: float,
        min: float = -np.inf,
        max: float = np.inf,
        vary: bool = True,
        flags: Tuple = (),
    ) -> None:
        """Add a parameter."""
        self.parameters[name] = FitParameter(name, value, min, max, not vary, flags)
        return self

    def set_value(self, name: str, value: float) -> None:
        """Set parameter value."""
        if name in self.parameters:
            self.parameters[name].value = value
        else:
            raise KeyError(f"Parameter {name} not found")

    def get_value(self, name: str) -> float:
        """Get parameter value."""
        return self.parameters[name].value

    def set_min(self, name: str, min: float) -> None:
        """Set minimum parameter value."""
        if name in self.parameters:
            self.parameters[name].min = min
        else:
            raise KeyError(f"Parameter {name} not found")

    def set_max(self, name: str, max: float) -> None:
        """Set maximum parameter value."""
        if name in self.parameters:
            self.parameters[name].max = max
        else:
            raise KeyError(f"Parameter {name} not found")

    def get_bounds(self) -> List[Tuple[float, float]]:
        """Get bounds for free parameters as list of (min, max) tuples."""
        return [(p.min, p.max) for p in self.parameters.values() if not p.fixed]

    def fix(self, name: str) -> None:
        """Fix parameter (don't vary during fit)."""
        self.parameters[name].fixed = True

    def free(self, name: str) -> None:
        """Free parameter (vary during fit)."""
        self.parameters[name].fixed = False

    def get_free_names(self) -> List[str]:
        """Get names of free (non-fixed) parameters."""
        return [name for name, p in self.parameters.items() if not p.fixed]

    def get_free_values(self) -> np.ndarray:
        """Get values of free parameters as array for scipy."""
        return np.array([p.value for p in self.parameters.values() if not p.fixed])

    def update_from_array(self, values: np.ndarray) -> None:
        """Update free parameter values from array (from scipy)."""
        free_params = [p for p in self.parameters.values() if not p.fixed]
        if len(values) != len(free_params):
            raise ValueError(
                f"Value array length {len(values)} doesn't match "
                f"free parameter count {len(free_params)}"
            )

        for param, value in zip(free_params, values):
            param.value = value

    def to_dict(self) -> Dict[str, float]:
        """Convert to dictionary {name: value}."""
        return {name: float(p.value) for name, p in self.parameters.items()}
    valuesdict = to_dict

    def copy(self):
        result = Parameters()
        for name, p in self.parameters.items():
            result.add(
                name,
                value=p.value,
                min=p.min,
                max=p.max,
                vary=not p.fixed,
                flags=p.flags,
            )
        return result

    def __repr__(self) -> str:
        lines = ["Parameters:"]
        for p in self.parameters.values():
            lines.append(f"  {p}")
        return "\n".join(lines)


class LMFitter:
    """
    Levenberg-Marquardt non-linear least squares fitter.
    """

    def __init__(
        self,
        ftol: float = 1e-9,
        xtol: float = 1e-9,
        gtol: float = 1e-9,
        max_nfev: int = 500,
        verbose: int = 0,
    ):
        """
        Initialize LM fitter.

        Parameters
        ----------
        ftol : float
            Tolerance for termination by change in cost function
        xtol : float
            Tolerance for termination by change in parameters
        gtol : float
            Tolerance for termination by gradient norm
        max_nfev : int
            Maximum number of function evaluations
        verbose : int
            Verbosity level: 0 (silent), values greater than 0 print a summary
        """
        self.ftol = ftol
        self.xtol = xtol
        self.gtol = gtol
        self.max_nfev = max_nfev
        self.verbose = verbose

    @staticmethod
    def _to_internal(values: np.ndarray, bounds: List[Tuple[float, float]]) -> np.ndarray:
        """Map physical parameter values to unconstrained LM variables.

        These are the Minuit-style transformations used by lmfit.  They let
        the Levenberg-Marquardt algorithm operate without explicit bounds.
        """
        internal = np.empty_like(values, dtype=float)

        for i, (value, (lower, upper)) in enumerate(zip(values, bounds)):
            if lower > upper:
                raise ValueError(f"Invalid bounds [{lower}, {upper}]")
            if value < lower or value > upper:
                raise ValueError(
                    f"Initial value {value} is outside bounds [{lower}, {upper}]"
                )

            if np.isneginf(lower) and np.isposinf(upper):
                internal[i] = value
            elif np.isposinf(upper):
                internal[i] = np.sqrt((value - lower + 1.0) ** 2 - 1.0)
            elif np.isneginf(lower):
                internal[i] = np.sqrt((upper - value + 1.0) ** 2 - 1.0)
            elif lower == upper:
                raise ValueError(
                    "A varying parameter cannot have identical lower and upper bounds"
                )
            else:
                scaled = 2.0 * (value - lower) / (upper - lower) - 1.0
                internal[i] = np.arcsin(np.clip(scaled, -1.0, 1.0))

        return internal

    @staticmethod
    def _from_internal(
        internal: np.ndarray, bounds: List[Tuple[float, float]]
    ) -> np.ndarray:
        """Map unconstrained LM variables back to physical parameter values."""
        values = np.empty_like(internal, dtype=float)

        for i, (value, (lower, upper)) in enumerate(zip(internal, bounds)):
            if np.isneginf(lower) and np.isposinf(upper):
                values[i] = value
            elif np.isposinf(upper):
                values[i] = lower - 1.0 + np.sqrt(value * value + 1.0)
            elif np.isneginf(lower):
                values[i] = upper + 1.0 - np.sqrt(value * value + 1.0)
            else:
                values[i] = lower + 0.5 * (np.sin(value) + 1.0) * (upper - lower)

        return values

    @staticmethod
    def _external_derivative(
        internal: np.ndarray, bounds: List[Tuple[float, float]]
    ) -> np.ndarray:
        """Return dx/du for conversion of covariance to physical space."""
        derivative = np.empty_like(internal, dtype=float)

        for i, (value, (lower, upper)) in enumerate(zip(internal, bounds)):
            if np.isneginf(lower) and np.isposinf(upper):
                derivative[i] = 1.0
            elif np.isposinf(upper):
                derivative[i] = value / np.sqrt(value * value + 1.0)
            elif np.isneginf(lower):
                derivative[i] = -value / np.sqrt(value * value + 1.0)
            else:
                derivative[i] = 0.5 * (upper - lower) * np.cos(value)

        return derivative

    def fit(
        self,
        residual_func: Callable,
        params: Parameters,
        args: Optional[Tuple] = (),
        num_data_points: int = -1,
    ) -> Tuple[Parameters, Dict]:
        """
        Perform Levenberg-Marquardt fit.

        Parameters
        ----------
        residual_func : callable
            Function that takes parameter array and returns residuals
            Signature: residual_func(param_array) -> np.ndarray
        params : Parameters
            Initial parameters (will be modified in place)
        args : tuple
            Optional arguments to residual_func

        Returns
        -------
        params : Parameters
            Fitted parameters
        result_info : dict
            Fit statistics and convergence information
        """
        bounds_list = params.get_bounds()

        # The optimizer sees unconstrained internal variables.  The residual
        # function and Parameters object continue to see physical values.
        def objective(internal, *args):
            values = self._from_internal(internal, bounds_list)
            params.update_from_array(values)
            return residual_func(params, *args)
        
        # Get initial values and bounds
        initial_values = params.get_free_values()
        if len(initial_values) > 0:
            x0 = self._to_internal(initial_values, bounds_list)

            # Bounds are enforced by the parameter transformations, so the
            # MINPACK Levenberg-Marquardt implementation can be used.  Its
            # internal diagonal scaling is important for parameters spanning
            # many orders of magnitude.
            best, cov_internal, info, message, status = leastsq(
                lambda internal: objective(internal, *args),
                x0,
                full_output=True,
                ftol=self.ftol,
                xtol=self.xtol,
                gtol=self.gtol,
                maxfev=self.max_nfev,
                epsfcn=1e-10,
                factor=100,
            )

            result = OptimizeResult(
                x=best,
                fun=info['fvec'],
                success=status in (1, 2, 3, 4),
                message=message,
                status=status,
                nfev=info['nfev'],
                cost=0.5 * np.sum(info['fvec'] ** 2),
                # qtf is MINPACK's available first-order convergence measure.
                optimality=np.linalg.norm(info.get('qtf', np.array([np.nan])), ord=np.inf),
                jac=None,
            )

            # Update parameters with final physical values.
            params.update_from_array(self._from_internal(result.x, bounds_list))

            # Compute uncertainties from Jacobian
            uncertainties = {}
            if num_data_points > len(result.x) and cov_internal is not None:
                try:
                    # leastsq returns the unscaled covariance in internal
                    # coordinates, equivalent to inv(J.T @ J).
                    pcov = cov_internal.copy()

                    # Compute the residual variance (Reduced Chi-Squared)
                    num_parameters = len(result.x)
                    degrees_of_freedom = num_data_points - num_parameters

                    # res.fun contains the residual array at the optimal solution
                    residual_sum_of_squares = np.sum(result.fun**2)
                    reduced_chi_sq = residual_sum_of_squares / degrees_of_freedom

                    # Final Covariance Matrix
                    pcov = pcov * reduced_chi_sq

                    # Convert internal-coordinate covariance to the physical
                    # parameter space.
                    dx_du = self._external_derivative(result.x, bounds_list)
                    pcov = pcov * np.outer(dx_du, dx_du)

                    # Standard errors = sqrt(diagonal of covariance)
                    std_errors = np.sqrt(np.diag(pcov))

                    free_names = params.get_free_names()
                    for name, std_err in zip(free_names, std_errors):
                        uncertainties[name] = std_err

                except np.linalg.LinAlgError:
                    # Singular matrix - uncertainties not computable
                    pass

        else:
            print("WARNING: zero model parameters to be optimized!")
            result = OptimizeResult()
            result.x = initial_values
            result.success = True
            result.message = 'No parameters to be fitted'
            result.nfev = 0
            result.cost = 0.0
            result.optimality = 0.0
            uncertainties = {}

        # Compute WSSR
        final_residuals = objective(result.x, *args)
        wssr = np.sum(final_residuals**2)
        rmse = np.sqrt(np.mean(final_residuals**2))

        # Package results
        result_info = {
            'success': result.success,
            'message': result.message,
            'nfev': result.nfev,
            'wssr': wssr,
            'rmse': rmse,
            'cost': result.cost,
            'optimality': result.optimality,
            'uncertainties': uncertainties,
        }

        if self.verbose > 0:
            print(f"\n=== Fit Results ===")
            print(f"Status: {result.message}")
            print(f"Iterations: {result.nfev}")
            print(f"WSSR: {wssr:.6e}")
            print(f"\nFitted parameters:")
            for name in params.get_free_names():
                value = params.get_value(name)
                if name in uncertainties:
                    unc = uncertainties[name]
                    rel_unc = (unc / value * 100) if value != 0 else 0
                    print(f"  {name} = {value:.6g} ± {unc:.4g} ({rel_unc:.2f}%)")
                else:
                    print(f"  {name} = {value:.6g}")

        return params, result_info


def minimize(
    residual_func: Callable,
    initial_params: Parameters,
    args: Optional[Tuple] = (),
    ftol: float = 1e-9,
    max_nfev: int = 500,
    verbose: int = 0,
) -> Tuple[Dict[str, float], Dict[str, float], Dict]:
    """
    Convenience function for fitting.

    Parameters
    ----------
    residual_func : callable
        Function to compute residuals: func(params_dict) -> residuals
    initial_params : Parameters
        Initial parameter values and bounds
    args : tuple
        Optional arguments to residual_func
    ftol : float
        Convergence tolerance
    max_nfev : int
        Maximum iterations
    verbose : int
        Verbosity level

    Returns
    -------
    fitted_params : dict
        Fitted parameter values
    uncertainties : dict
        Parameter uncertainties
    fit_info : dict
        Additional fit information
    """
    fitter = LMFitter(ftol=ftol, max_nfev=max_nfev, verbose=verbose)
    fitted_params, result_info = fitter.fit(residual_func, initial_params, args=args)

    result_info['params'] = fitted_params.to_dict()
    return result_info
