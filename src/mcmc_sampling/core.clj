(ns mcmc-sampling.core
  (:require [mcmc-sampling.util :as util]
            [mcmc-sampling.gibbs :as gibbs]
            [clojure.core.matrix :as m])
  (:gen-class))

(def mu1 [0 0])
(def mu2 [1 1])
(def a1 [[1 0] [0 1]])
(def a2 [[1 0] [0 1]])
(def number-of-samples 100000)
(def d 0.1)
(def area-of-d (* Math/PI (util/sqr d)))

(def f1 (partial util/gauss-function-2d mu1 a1))
(def f2 (partial util/gauss-function-2d mu2 a2))

(defn multi-gauss-function
  "Returns the average of two gauss functions. "
  [x]
  (/ (+ (f1 x) (f2 x)) 2))

(def density-function multi-gauss-function)

(defn- mh-get-next-sample
  "Get next sample in MH algorithm"
  [previous-sample]
  (let [new-sample (util/gauss-sample-2d)
        threshold (/ (* (util/standard-gauss-function-2d previous-sample)
                        (density-function new-sample))
                     (* (util/standard-gauss-function-2d new-sample)
                        (density-function previous-sample)))
        alpha (min threshold 1)
        u (util/uniform-sample)]
    (if (<= u alpha)
      new-sample
      previous-sample)))

(defn- mh-get-samples
  "Get samples using MH algorithm."
  [f]
  (loop [samples [[0 0]]
         i 0]
    (if (= i number-of-samples)
      samples
      (recur (conj samples
                   (mh-get-next-sample (peek samples)))
             (inc i)))))

(defn- mc-calc-probability-density
  "Calc probability density of point x using Monte Carlo Method."
  [x]
  (let [samples (mh-get-samples density-function)]
    (loop [i 0
           cnt 0]
      (if (= i number-of-samples)
        (/ (/ cnt number-of-samples)
           area-of-d)
        (recur (inc i)
               (+ cnt (util/neighbor? x (samples i) d)))))))


(def run mc-calc-probability-density)

(defn -main
  []
  (println (run [0.1 0.1])))
