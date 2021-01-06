# fft-kotlin
dynamic time warping - in kotlin


            val distFn = DistanceFunctionFactory.getDistFnByName("EuclideanDistance");
            var info = DTW.getWarpInfoBetween(TimeSeries(Rarr2), TimeSeries(Rarr8), distFn);
            
            info.distance
            info.path
