{-# LANGUAGE 
        StandaloneDeriving,
        FlexibleContexts,
        UndecidableInstances,
        RecordWildCards
  #-}
module Math.Quadrature.Gaussian
    ( Estimate(..)
    , QRule(..)
    , integrate, integrateRange
    , gaussLegendre
    ) where

import Math.Quadrature.Simple (Estimate(..))

import qualified Data.Vector as V
import qualified Data.Vector.Generic as GV
import qualified Data.Vector.Generic.Mutable as MV
import Math.Polynomial
import Math.Polynomial.Legendre

data QRule v a = QRule
    { qRange    :: !(a,a)
    , qTable    :: !(v (a,a))
    }

deriving instance (Eq   a, Eq   (v (a,a))) => Eq   (QRule v a)
deriving instance (Show a, Show (v (a,a))) => Show (QRule v a)

-- |Apply a Gaussian quadrature rule to a function.  The endpoints of integration
-- are a property of the rule - typically, they were supplied to the function
-- (such as 'gaussLegendre' below) that computed the parameters for the rule.
integrate :: (GV.Vector v (a,a), Num a) => QRule v a -> (a -> a) -> Estimate a
integrate QRule{..} f = Estimate n Nothing $ sum [w * f x | (x, w) <- GV.toList qTable]
    where n = toInteger (GV.length qTable)

-- |Apply a Gaussian quadrature rule to a function, using user-specified 
-- endpoints of integration which may differ from those for which the
-- quadrature rule was originally defined.
integrateRange :: (GV.Vector v (a,a), Fractional a) => QRule v a -> (a -> a) -> a -> a -> Estimate a
integrateRange qRule f a0 a1
    | x0 == a0 && x1 == a1  = integrate qRule f
    | otherwise             = integrate qRule f'
    where
        (x0, x1) = qRange qRule
        scale  = (a1 - a0) / (x1 - x0)
        
        f' x = scale * f ((x - x0) * scale + a0)

-- |Given the lower and upper limits of integration, number of points n, and
-- desired accuracy eps, this routine returns a Gaussian quadrature rule containing
-- the abscissas and weights of the Gauss-Legendre n-point quadrature formula.
--
-- Typical values for eps are 1e-15 for Double, 1e-6 for Float.  The specified
-- eps will be approximately the accuracy of the integration method as well,
-- for suitable integrands.
gaussLegendre :: (GV.Vector v a, GV.Vector v (a,a), Fractional a, Ord a)
    => a -> a -> Int -> a -> QRule v a
gaussLegendre x1 x2 n eps = QRule (x1, x2) table
    where
        roots = GV.fromList (legendreRoots n eps)
        derivs = GV.map (snd . evalLegendreDeriv n) roots
        
        table = flip GV.map (GV.zip roots derivs) $ \(z,dy) ->
            ( {- abscissa -}    xm + xl * z
            , {- weight   -}    2 * xl / ((1 - z*z)*dy*dy)
            ) where
                xm = 0.5*(x2+x1)
                xl = 0.5*(x2-x1)
