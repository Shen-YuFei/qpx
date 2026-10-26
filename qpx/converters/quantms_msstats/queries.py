"""SQL for normalized QuantMS MSstats rows."""

_NORMALIZED_SQL = """
            CREATE OR REPLACE TEMP TABLE _normalized AS
            SELECT trim(CAST(m.$protein AS VARCHAR)) AS protein_raw,
                   trim(CAST(m.$peptide AS VARCHAR)) AS peptide_raw,
                   p.sequence,
                   p.peptidoform,
                   TRY_CAST(m.$charge AS SMALLINT) AS charge,
                   COALESCE(rr.run_file_name, ru.run_file_name) AS run_file_name,
                   rr.run_file_name AS reference_file,
                   ru.run_file_name AS run_file,
                   TRY_CAST(m.$intensity AS FLOAT) AS intensity,
                   $rt_expr AS rt,
                   $scan_expr AS scan_number,
                   $reference_expr AS reference_raw,
                   $label_expr AS intensity_label,
                   $order_expr AS channel_order,
                   count(DISTINCT trim(CAST(m.$protein AS VARCHAR)))
                     OVER (PARTITION BY p.sequence) = 1
                     AND NOT contains(trim(CAST(m.$protein AS VARCHAR)), ';') AS is_unique
            FROM _msstats_input m
            LEFT JOIN _run_aliases rr ON rr.alias_key = $reference_key
            LEFT JOIN _run_aliases ru ON ru.alias_key = $run_key
            LEFT JOIN _peptide_lookup p
              ON p.peptide_raw = trim(CAST(m.$peptide AS VARCHAR))
            $channel_join
            """
