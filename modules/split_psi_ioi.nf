process split_psi_ioi {
    publishDir "${params.outdir}/suppa/split_psi_ioi", mode: 'copy'

    input:
    path(psi_file)
    path(samplesheet)

    output:
    path "split_ioi_*.psi", emit: split_psi_ioi

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import os

    # Leer el samplesheet
    samples = pd.read_csv("${samplesheet}", sep=',')
    sample_to_condition = dict(zip(samples['sample_id'], samples['condition']))

    # Leer el archivo PSI sin interpretar los NA/NaN
    psi = pd.read_csv("${psi_file}", sep='\t', index_col=0, na_filter=False)

    # Verificar que todas las muestras del samplesheet estén en el archivo PSI
    missing_samples = set(sample_to_condition.keys()) - set(psi.columns)
    if missing_samples:
        print(f"Advertencia: Muestras faltantes en el archivo PSI IOI: {missing_samples}")

    # Separar por condición conservando la estructura del archivo original
    for condition in set(sample_to_condition.values()):
        samples_in_condition = [s for s, c in sample_to_condition.items() if c == condition and s in psi.columns]
        if samples_in_condition:
            output_file = f"split_ioi_{os.path.basename('${psi_file}').replace('.psi', '')}_{condition}.psi"
            # Crear un DataFrame con la estructura original y valores para las muestras de esta condición
            split_psi = psi[samples_in_condition].copy()
            # Guardar el archivo manteniendo los valores NA/NaN exactamente como están
            with open(output_file, 'w') as f:
                # Escribir el encabezado con los nombres de las muestras
                f.write("\t".join(samples_in_condition) + "\\n")
                # Escribir los datos línea por línea para mantener el formato exacto
                for index, row in split_psi.iterrows():
                    values = "\t".join(row.astype(str))
                    f.write(f"{index}\t{values}\\n")
            print(f"Archivo IOI generado: {output_file}")
        else:
            print(f"Advertencia: No se encontraron muestras para la condición {condition} en el archivo PSI IOI.")
        """
    }
