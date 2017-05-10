(ns mcmc-sampling.core
  (:require [mcmc-sampling.util :as util]
            [clojure.core.matrix :as m])
  (:gen-class))

(def mu1 [0 0])
(def a1 [[1 0] [0 1]])
(def mu2 [1 1])
(def a2 [[1 0] [0 1]])
(def number-of-samples 100000)
(def d 0.1)
(def area-of-d (* Math/PI (util/sqr d)))

(def f1 (partial util/gauss-function-2d mu1 a1))
(def f2 (partial util/gauss-function-2d mu2 a2))
;; (def f2 f1)

(defn density-function
  [x]
  (/ (+ (f1 x) (f2 x)) 2))

(defn- neighbor?
  [x y]
  (let [dis (util/distance x y)]
    (if (< dis d) 1 0)))

(def get-sample identity)

(defn- gibbs-get-sample
  []
  (let [partial-x (fn [x] #(density-function [x %]))
        partial-y (fn [y] #(density-function [% y]))
        burn-in 20]
    (loop [i 0
           x 0
           y 0]
      (if (= i burn-in)
        [x y]
        (let [x (get-sample (partial-x x))
              y (get-sample (partial-y y))]
          (recur (inc i) x y))))))

(defn- MH-get-next-sample
  [previous-sample]
   (let [new-sample (util/gauss-sample-2d)
        threshold (/ (* (standard-gauss-function previous-sample)
                        (density-function new-sample))
                     (* (standard-gauss-function new-sample)
                        (density-function previous-sample)))
        alpha (min threshold 1)
        u (util/uniform-sample)]
    (if (<= u alpha)
      new-sample
      previous-sample)))

(defn- MH-get-samples
  [f]
  (loop [samples [[0 0]]
         i 0]
    (if (= i number-of-samples)
      samples
      (recur (conj samples
                   (MH-get-next-sample (peek samples)))
             (inc i)))))

(defn- mc-calc-probability-density
  [x]
  (let [samples (MH-get-samples density-function)]
    (loop [i 0
           cnt 0]
      (if (= i number-of-samples)
        (/ (/ cnt number-of-samples)
           area-of-d)
        (recur (inc i)
               (+ cnt (neighbor? x (samples i))))))))


(def run mc-calc-probability-density)

(def -main #())
