SIMULACION Y MODELOS NULOS

1. Simulación. Configurar el entorno de simulacion con el fichero config_values.R
2. Simulación. Lanzar el simulador null_models_spectral_analysis que genera los directorios de resultados en XXXsmodels siendo XXX el valor de debugoref. Si el parámetro cold_start es TRUE se borran todos los resultados anteriores. El programa pregunta antes de borrar
3. Análisis. Lanzar normalize_magnitudes que crea los ficheros ALLNORMALIZED.csv y networktypes.csv
4. Análisis. Lanzar compute_stats_normalized_distances.R que crea los ficheros de correlaciones y el directorio "normnalized" dentro de "results"
5. Gráficas. plot_normdistavgs.R, genera el directorio "analysis" dentro de "plots" 
5. Gráficas. plot_matrix.R, genera el directorio "analysis" dentro de "plots" 


ANALISIS DE ANIDAMIENTO

1. Simulación, análisis y gráficas eigen_analysis.R  No puede lanzarse si no se ha completado una simulación completa con null_models_spectral_analysis porque lee los datos de "NetworkMagnitudes.csv"

