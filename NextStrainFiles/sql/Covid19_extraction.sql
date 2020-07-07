with v as (select NOLSPQ as no_lspq, AGE_ANNEE as age,SEX as sex, RSS_PATIENT as rss_patient,AUCUN_VOYAGE as aucun_voyage, VOYAGE_PAYS_1 as voyage_pays_1, DATE_PRELEVEMENT as date_prelev, DATE_RECEPTION as date_recu,
CH as ch,RESULTAT_LABORATOIRE as res_lab, SUBSTR(POSTAL_CODE,1,3) as postal_code ,MAX(nCov2019_nCov_ARN) as max_res, 

case
    when Max(nCoV2019_nCoV_Ct) is null
    then '0'
    else TO_CHAR(ROUND(Max(nCoV2019_nCoV_Ct)))
end as max_ct,

ADDR_CH as addr_ch

from (select distinct f.folderno as NOLSPQ, to_char(cr.BIRTH_DATE,'YYYY-MM-DD') as DATE_NAISSANCE,
p.LAST_NAME as NOM, p.FIRST_NAME as PRENOM, 
COALESCE(case 
    when FLOOR(MONTHS_BETWEEN(COALESCE(cr.DATE_COLLECTED,cr.DATE_RECEIVED),cr.BIRTH_DATE )/12) > 1
    THEN TRUNC(MONTHS_BETWEEN(COALESCE(cr.DATE_COLLECTED,cr.DATE_RECEIVED),cr.BIRTH_DATE)/12,0)
    else
        TRUNC(MONTHS_BETWEEN(COALESCE(cr.DATE_COLLECTED,cr.DATE_RECEIVED),cr.BIRTH_DATE )/12,2)
end,999) as AGE_ANNEE,
COALESCE(p.SEX,'NA') as SEX, 
 
case
    when cr.OWNER_COUNTY is null
    then (select REGION_NAME from REGION where TYPE = 'RSS' and REGION_CODE = rc.LSPQ_RSS)
    else cr.OWNER_COUNTY
end as RSS_PATIENT,

COALESCE(cr.OWNER_ZIP,rc.ZIP) as POSTAL_CODE,
case 
    when crm.FIELD45 is null 
    then 'O' 
    else 'N' 
end AS AUCUN_VOYAGE,
case 
    when crm.FIELD45 is null
    then 'AUCUN_VOYAGE'
    else LOWER(crm.FIELD45)
end as VOYAGE_PAYS_1,

nvl(to_char(cr.DATE_COLLECTED, 'YYYY-MM-DD'),'') AS DATE_PRELEVEMENT,
nvl(to_char(cr.DATE_RECEIVED, 'YYYY-MM-DD'),'') AS DATE_RECEPTION,
rc.COMPNAME AS CH,

(
select REPLACE(REPLACE(max(concat(concat(concat(concat(concat(adress,', '),city),', '),state),', Candada')),CHR(10),' '),CHR(13),' ') FROM RASCLIENTS where compname = rc.COMPNAME
group by compname)as ADDR_CH,

case 
   when replace(replace(replace(crm.FIELD06, chr(13), ' '), chr(10),' '),';',':') is null
   then ' '
   else replace(replace(replace(crm.FIELD06, chr(13), ' '), chr(10),' '),';',':')
end as RESULTAT_LABORATOIRE,

( select    
    Final
    from RESULTS coro
    where coro.ORDNO = ot.ORDNO  AND coro.TESTNO='2019-nCoV'  AND coro.ANALYTE='2019-nCoV ARN'  and coro.REPORTABLE = 'Y' 
    and coro.ORIGREC = ( select max(origrec) from results res2 where res2.ordno=coro.ordno and res2.testno=coro.testno AND res2.TESTNO='2019-nCoV' and res2.ANALYTE='2019-nCoV ARN' )
) as nCoV2019_nCoV_ARN,

( select 
    case 
        when res.Final > 100
        then '0'
        else res.Final 
    end
    
    from RESULTS res
    where res.ORDNO = ot.ORDNO  AND res.TESTNO='2019-nCoV'  AND res.ANALYTE='CoVN (LSPQ) Ct'  
    and res.ORIGREC = ( select max(origrec) from results res2 where res2.ordno=res.ordno and res2.testno=res.testno AND res2.TESTNO='2019-nCoV' and res2.ANALYTE='CoVN (LSPQ) Ct' )
) as nCoV2019_nCoV_Ct

from 
		CENTRALRECEIVING cr 
		inner join FOLDERS f on f.FOLDERNO = cr.EXTERNAL_ID 
		inner join METADATA crm on crm.ID = cr.METADATA_GUID 
		inner join RASCLIENTS rc on rc.RASCLIENTID = cr.RASCLIENTID 
		left join PATIENTS p on p.PID = cr.PID 
		inner join ORDERS o on o.FOLDERNO = cr.EXTERNAL_ID 
		inner join ORDTASK ot on ot.ORDNO = o.ORDNO
        
where 
		ot.TESTCODE in (2685,2689,2167) and cr.DATE_RECEIVED > to_date('2020-01-01', 'YYYY-MM-DD') 
        --and cr.PANEL_LIST like '2019-nCoV%' and rc.RASCLIENTID not in ('LSPQCEC','LSPQCIC','LSPQF','LSPQP','LSPQV') 
        and (cr.PANEL_LIST like '2019-nCoV%' or cr.PANEL_LIST like '%COVID-19 - Confirmation%' or cr.PANEL_LIST like '%COVID-19 - Soustraitant%') and rc.RASCLIENTID not in ('LSPQCEC','LSPQCIC','LSPQF','LSPQP','LSPQV') 
        
order by f.FOLDERNO)   GROUP BY NOLSPQ, AGE_ANNEE,SEX, RSS_PATIENT,AUCUN_VOYAGE, VOYAGE_PAYS_1, DATE_PRELEVEMENT, DATE_RECEPTION,
CH,RESULTAT_LABORATOIRE, SUBSTR(POSTAL_CODE,1,3),ADDR_CH ORDER BY NOLSPQ) select v.no_lspq,v.age,v.sex,v.rss_patient,v.aucun_voyage, v.voyage_pays_1, v.date_prelev,
v.date_recu, v.ch, v.postal_code,  v.max_res, v.max_ct,v.res_lab,v.addr_ch from v where  v.max_res in ('Détecté','détecté'); 




