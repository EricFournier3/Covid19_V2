select distinct FOLDERNO as NUMERO_SGIL , to_char(cr.DATE_COLLECTED,'YYYY-MM-DD') as SAMPLED_DATE, cr.LAST_NAME as NOM, cr.FIRST_NAME as PRENOM, to_char(cr.BIRTH_DATE,'YYYY-MM-DD') as DATE_NAISS, PID 
from CENTRALRECEIVING cr inner join RASCLIENTS rc on rc.RASCLIENTID = cr.RASCLIENTID where (cr.PANEL_LIST like '2019-nCoV%' or cr.PANEL_LIST like '%COVID-19 - Confirmation%' or  cr.PANEL_LIST like '%GeneXpert SARS-CoV-2%' or
cr.PANEL_LIST like '%COVID-19 - Soustraitant%' or cr.PANEL_LIST like '%COVID-19 - Surveillance WGS%') and rc.RASCLIENTID
not in ('LSPQCEC','LSPQCIC','LSPQF','LSPQP','LSPQV') AND FOLDERNO  in (SELECT FOLDERNO from RESULTS res where (res.TESTCODE in (2705) and res.ANALYTE='Génotype') or (res.TESTCODE in (2685,2689)
and res.ANALYTE='2019-nCoV ARN' and res.Final in ('Détecté','détecté')) or (res.TESTCODE in (2704)) and res.ANALYTE='SARS-CoV-2' and Final = 'Détecté' );

