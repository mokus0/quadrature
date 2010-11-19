{-# LANGUAGE BangPatterns #-}
-- |Simple quadrature schemes based on repeated subdivision, implemented
-- as functions returning series of estimates of the integral.
-- 
-- There is much room for stylistic and usability improvements.  Most notably,
-- some convenient drivers should be written and integrand transformations should
-- be implemented.  Also, it would be nice for the domain of integration to be
-- generalized to something like range-sets.
module Math.Quadrature.Simple where

import Data.Bits
import Data.List
import Data.Maybe (fromMaybe)
import Math.Polynomial.Interpolation

data Estimate a = Estimate
    { numPoints :: Integer
    , stepSize  :: Maybe a  -- 'Maybe' because not all algorithms use equally spaced points
    , estimate  :: a
    } deriving (Eq, Ord, Show)

-- |@trapezoidRule f a b@ computes a sequence of improving estimates of the
-- integral of f from a to b, using the trapezoid rule.
trapezoidRule :: (Fractional b, Ord b) => (b -> b) -> b -> b -> [Estimate b]
trapezoidRule f a b = iterate next first
    where
        first = let h0 = (b - a)
                 in Estimate 2 (Just h0) (0.5 * h0 * (f a + f b))
        next (Estimate !nPts (Just !h0) i0) = Estimate (nPts+nPts-1) (Just h1) i1
            where
                h1 = 0.5 * h0
                i1 = 0.5 * i0 + h1 * foldl1' (+) ys
                ys = map f (midpoints h0)
        
        -- |Divide the interval [a,b] into segments of width dx, and return
        -- the midpoint of each segment
        midpoints dx = takeWhile (< max a b) (iterate (+dx) x0)
            where
                x0 = min a b + 0.5 * abs dx

-- |@simpson'sRule f a b@ computes a sequence of improving estimates of the
-- integral of f from a to b, using Simpson's rule.
simpson'sRule :: (Fractional b, Ord b) => (b -> b) -> b -> b -> [Estimate b]
simpson'sRule f a b = 
    [ Estimate n h ((4/3) * s_2n - (1/3) * s_n)
    | Estimate _ _ s_n : Estimate n h s_2n : _ <- tails (trapezoidRule f a b)
    ]


midpointRule :: (Fractional a) => (a -> a) -> a -> a -> [Estimate a]
midpointRule f a b = iterate next first
    where
        h0 = b-a
        first = Estimate 1 (Just h0) (h0 * f (0.5 * (a+b)))
        extension x0 hop skip = x0 : extension (x0 + skip) skip hop
        next (Estimate n0 _ i) = Estimate n1 (Just h) (i / 3 + h * sum ys)
            where
                n1 = 3 * n0
                h  = h0 / fromInteger n1
                x0 = a + 0.5 * h
                ys = take (2 * fromInteger n0) (map f (extension x0 h (h+h)))

-- |@romberg k qrule f a b@ uses Romberg extrapolation to improve the 
-- series of estimates made by a quadrature rule by fitting an order-@k@ 
-- polynomial to (stepsize^2, X) and extrapolating to stepsize=0.  If the
-- quadrature rule provided does not populate the stepsize field, it is
-- estimated as the size of the interval divided by the number of points
-- sampled.
romberg :: (Fractional a) => Int -> (t -> a -> a -> [Estimate a]) -> t -> a -> a -> [Estimate a]
romberg k qrule f a b = 
    [ Estimate n h (polyInterp (take k terms) 0)
    | estimates <- tails (qrule f a b)
    , let terms = take k [(fromMaybe (guessH n) mbH ^ 2, x) | Estimate n mbH x <- estimates]
          Estimate n h _ = last (take k estimates)
    ]
    where
        guessH n = (b-a) / fromIntegral n

