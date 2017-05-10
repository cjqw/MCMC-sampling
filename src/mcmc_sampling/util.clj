(ns mcmc-sampling.util
  (:require [clojure.core.matrix :as m])
  (:gen-class))

(def sqr (fn [x] (* x x)))

(def uniform-sample rand)

(defn gauss-sample
  "Return a variable which satisfying gauss distribution."
  []
  (let [u (uniform-sample)
        v (uniform-sample)]
    (* (Math/sqrt (* -2 (Math/log u)))
       (Math/cos (* 2 Math/PI v)))))

(defn gauss-function
  [mu sigma x]
  (let [left (/ 1 (Math/sqrt (* 2 Math/PI)))
        index (/ (sqr (- x mu))
                 (* -2 (sqr sigma)))
        right (Math/exp index)]
    (* left right)))

(defn determinant
  "Calculate the determinant of a 2*2 matrix."
  [mat]
  (let [a (mat 0)
        b (mat 1)]
    (- (* (a 0) (b 1))
       (* (a 1) (b 0)))))

(defn gauss-function-2d
  "sigma = (A^{T} A)^{-1}"
  [mu a x]
  (let [v (m/mmul a (m/sub x mu))
        left (/ 1 (* 2 Math/PI))
        index (* -0.5 (m/mmul v v))
        right (Math/exp index)]
    (* left right)))

(defn gauss-sample-2d
  "Return a vector satisfying gauss-distribution
  where mu = [0 0] and a = [[1 0] [0 1]]."
  []
  (let [x (gauss-sample)
        y (gauss-sample)]
    [x y]))

(defn distance
  "Calculate the distance of two points."
  [[x0 y0] [x1 y1]]
   (Math/sqrt (+ (sqr (- x0 x1))
                 (sqr (- y0 y1)))))
