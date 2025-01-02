process split_tpm_file {
    publishDir "${params.outdir}/suppa/split_tpm", mode: 'copy'

    input:
    path(suppa_tpm)
    path(samplesheet)

    output:
    path "tpm_*.txt", emit: split_tpm
    path "conditions.txt", emit: conditions

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import os

    # Leer el samplesheet
    samples = pd.read_csv("${samplesheet}", sep=',')
    sample_to_condition = dict(zip(samples['sample_id'], samples['condition']))

    # Leer el archivo TPM preservando los tipos de datos originales
    tpm = pd.read_csv("${suppa_tpm}", sep='\\t', index_col=0, dtype=str)

    # Agrupar por condición
    condition_groups = {}
    for sample, condition in sample_to_condition.items():
        if sample in tpm.columns:
            if condition not in condition_groups:
                condition_groups[condition] = []
            condition_groups[condition].append(sample)
        else:
            print(f"Advertencia: La muestra {sample} no está en el archivo TPM.")

    # Separar por condición conservando la estructura del archivo original
    for condition, samples in condition_groups.items():
        output_file = f"tpm_{condition}.txt"
        # Crear un DataFrame con la estructura original y valores para las muestras de esta condición
        split_tpm = tpm[samples].copy()

        # Guardar el archivo manteniendo el formato exacto
        with open(output_file, 'w') as f:
            # Escribir el encabezado con los nombres de las muestras
            f.write("\\t".join(samples) + "\\n")
            # Escribir los datos línea por línea para mantener el formato exacto
            for index, row in split_tpm.iterrows():
                values = "\\t".join(row.astype(str))
                f.write(f"{index}\\t{values}\\n")

    # Crear archivo de condiciones
    with open("conditions.txt", "w") as f:
        for condition in set(sample_to_condition.values()):
            f.write(f"{condition}\\n")
    """
}
