process prepare_tpm_for_suppa {
    cache false
    publishDir "${params.outdir}/suppa/prepared_tpm", mode: 'copy'
    
    input:
    path tpm
    
    output:
    path "prepared_tpm.txt"
    path "tpm_preparation.log"
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    
    # Log file
    with open("tpm_preparation.log", "w") as log:
        try:
            log.write("=== Preparando archivo TPM para SUPPA ===\\n")
            
            # Leer el archivo TPM con manejo específico de espacios
            df = pd.read_csv("${tpm}", sep='\\s+', index_col=0)
            
            log.write(f"Dimensiones originales: {df.shape}\\n")
            log.write(f"Columnas: {', '.join(df.columns)}\\n")
            
            # Verificar y limpiar el índice (transcript IDs)
            log.write(f"Número de índices únicos: {len(df.index.unique())}\\n")
            
            # Reemplazar valores no numéricos y NaN con 0
            df = df.replace('\\s+', '0', regex=True)
            df = df.fillna(0)
            
            # Convertir todo a números
            for col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
            
            # Asegurar que todos los valores son ≥ 0
            df[df < 0] = 0
            
            # Guardar con formato consistente
            with open("prepared_tpm.txt", 'w') as f:
                # Escribir el encabezado
                header = df.columns.tolist()
                f.write("\\t".join(['transcript_id'] + header) + '\\n')
                
                # Escribir los datos
                for idx, row in df.iterrows():
                    values = [f"{x:.6f}" for x in row.values]
                    f.write(f"{idx}\\t" + "\\t".join(values) + '\\n')
            
            log.write("\\nArchivo TPM preparado guardado exitosamente\\n")
            log.write(f"Dimensiones finales: {df.shape}\\n")
            
        except Exception as e:
            log.write(f"\\nError durante la preparación: {str(e)}\\n")
            raise e
    """
}
