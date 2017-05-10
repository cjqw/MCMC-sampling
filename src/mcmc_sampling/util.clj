(ns mcmc-sampling.util
  (:require [clojure.core.matrix :as m])
  (:gen-class))

(def sqr (fn [x] (* x x)))

(defn determinant-2d
  "Calculate the determinant of a 2*2 matrix."
  [mat]
  (let [a (mat 0)
        b (mat 1)]
    (- (* (a 0) (b 1))
       (* (a 1) (b 0)))))

(def uniform-sample rand)

(defn gauss-sample
  "Return a variable which satisfying gauss distribution."
  [mu sigma]
  (let [u (uniform-sample)
        v (uniform-sample)
        x (* (Math/sqrt (* -2 (Math/log u)))
           (Math/cos (* 2 Math/PI v)))]
    (+ mu (* x sigma))))

(def standard-gauss-sample (partial gauss-sample 0 1))

(defn gauss-function
  "The density function of normal distribution."
  [mu sigma x]
  (let [left (/ 1 (Math/sqrt (* 2 Math/PI)))
        index (/ (sqr (- x mu))
                 (* -2 (sqr sigma)))
        right (Math/exp index)]
    (* left right)))

(defn gauss-function-2d
  "sigma = (A^{T} A)^{-1}"
  [mu a x]
  (let [v (m/mmul a (m/sub x mu))
        left (/ 1 (* 2 Math/PI))
        index (* -0.5 (m/mmul v v))
        right (Math/exp index)]
    (* left right)))

(defn standard-gauss-function-2d
  "Standard gauss function in 2d.
  That is, mu = [0 0] and a = [[1 0] [0 1]]."
  [x]
  (* (/ 1 (* 2 Math/PI))
     (Math/exp (m/mmul -0.5 x x))))

(defn gauss-sample-2d
  "Return a vector satisfying gauss-distribution
  where mu = [0 0] and a = [[1 0] [0 1]]."
  []
  (let [x (standard-gauss-sample)
        y (standard-gauss-sample)]
    [x y]))

(defn distance
  "Calculate the distance of two points."
  [[x0 y0] [x1 y1]]
   (Math/sqrt (+ (sqr (- x0 x1))
                 (sqr (- y0 y1)))))

(defn neighbor?
  "If the distance between x and y is less than d,
  return 1.Otherwise return 0."
  [x y d]
  (let [dis (distance x y)]
    (if (< dis d) 1 0)))
