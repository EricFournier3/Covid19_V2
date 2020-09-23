--select FOLDERNO,FINAL,ANALYTE from results where testcode = 2685 and ANALYTE = '2019-nCoV ARN' order by folderno,analyte;
--select FOLDERNO,FINAL,ANALYTE from results where testcode = 2689 and ANALYTE = '2019-nCoV ARN' order by folderno,analyte;

--SELECT DISTINCT FINAL FROM RESULTS WHERE  testcode = 2685 and ANALYTE = '2019-nCoV ARN' ;
--SELECT DISTINCT FINAL FROM RESULTS WHERE  testcode = 2689 and ANALYTE = '2019-nCoV ARN' ;
--SELECT DISTINCT FINAL FROM RESULTS WHERE  testcode = 2704 and ANALYTE = 'SARS-CoV-2' ;

select NUMERO_SGIL,DATE_NAISS,NOM,PRENOM,AGE_ANNEE,SEX,NAM,PID,RSS_PATIENT,AUCUN_VOYAGE,TRAVEL_HISTORY,SAMPLED_DATE,CH_NAME,RSS_CH,
CASE WHEN nCoV_LNM2019_nCoV_ARN is NULL then 'NO RESULT' ELSE nCoV_LNM2019_nCoV_ARN END AS LNM_RESULT,
case when nCoV_LSPQ2019_nCoV_ARN is null then 'NO RESULT' ELSE nCoV_LSPQ2019_nCoV_ARN END AS LSPQ_RESULT,
case when geneXpertSARScov2 is null then 'NO RESULT' ELSE geneXpertSARScov2 END as GENEXPERT_RESULT,
case when WGS_COVID is null then 'NO RESULT' ELSE WGS_COVID END as  WGS_COVID_RESULT,
CT,
CH_ADRESS,
POSTAL_CODE

FROM (SELECT DISTINCT f.FOLDERNO as NUMERO_SGIL, to_char(cr.BIRTH_DATE, 'YYYY-MM-DD') as DATE_NAISS, 
		p.LAST_NAME AS NOM, 
		p.FIRST_NAME AS PRENOM,
        COALESCE(case 
          when FLOOR(MONTHS_BETWEEN(COALESCE(cr.DATE_COLLECTED,cr.DATE_RECEIVED),cr.BIRTH_DATE )/12) > 1
         THEN TRUNC(MONTHS_BETWEEN(COALESCE(cr.DATE_COLLECTED,cr.DATE_RECEIVED),cr.BIRTH_DATE)/12,0)
         else
        TRUNC(MONTHS_BETWEEN(COALESCE(cr.DATE_COLLECTED,cr.DATE_RECEIVED),cr.BIRTH_DATE )/12,2)
        end,999) as AGE_ANNEE,
    COALESCE(p.SEX,'NA') as SEX,
    p.SSN AS NAM,
    cr.PID AS  PID,
    
    case
      when cr.OWNER_COUNTY is null
      then (select REGION_NAME from REGION where TYPE = 'RSS' and REGION_CODE = rc.LSPQ_RSS)
      else cr.OWNER_COUNTY
   end as RSS_PATIENT,
   
   COALESCE(cr.OWNER_ZIP,rc.ZIP) as POSTAL_CODE,
    
    case when crm.FIELD45 is null then 'O' else 'N' end AS AUCUN_VOYAGE,  
    
    case 
      when crm.FIELD45 is null
      then 'AUCUN_VOYAGE'
      else LOWER(crm.FIELD45)
    end as TRAVEL_HISTORY,
    
    nvl(to_char(cr.DATE_COLLECTED, 'YYYY-MM-DD'),'') AS SAMPLED_DATE,
    rc.COMPNAME AS CH_NAME,
    rc.LSPQ_RSS AS RSS_CH,
    --LNM
    (SELECT coroLNM.FINAL
      FROM RESULTS coroLNM
       WHERE coroLNM.ORDNO = ot.ORDNO and  coroLNM.TESTNO = '2019-nCoV - LNM' and coroLNM.ANALYTE='2019-nCoV ARN' and coroLNM.REPORTABLE = 'Y'
       AND coroLNM.ORIGREC = ( select max(origrec) from results res2 where res2.ordno=coroLNM.ordno and res2.testno=coroLNM.testno AND res2.TESTNO='2019-nCoV - LNM' and res2.ANALYTE='2019-nCoV ARN')
     )AS nCoV_LNM2019_nCoV_ARN,
     
     --LSPQ
    (SELECT 
      coroLSPQ.FINAL
      FROM RESULTS coroLSPQ
       WHERE coroLSPQ.ORDNO = ot.ORDNO and  coroLSPQ.TESTNO = '2019-nCoV' and coroLSPQ.ANALYTE='2019-nCoV ARN' and coroLSPQ.REPORTABLE = 'Y'
       AND coroLSPQ.ORIGREC = ( select max(origrec) from results res2 where res2.ordno=coroLSPQ.ordno and res2.testno=coroLSPQ.testno AND res2.TESTNO='2019-nCoV' and res2.ANALYTE='2019-nCoV ARN')
     )AS nCoV_LSPQ2019_nCoV_ARN,
     
     --GENEXPERT
    (SELECT 
      geneXpert.FINAL
      FROM RESULTS geneXpert
       WHERE geneXpert.ORDNO = ot.ORDNO and  geneXpert.TESTNO = 'Xpert Xpress SARS-CoV-2' and geneXpert.ANALYTE='SARS-CoV-2' and geneXpert.REPORTABLE = 'Y'
       AND geneXpert.ORIGREC = ( select max(origrec) from results res2 where res2.ordno=geneXpert.ordno and res2.testno=geneXpert.testno AND res2.TESTNO='Xpert Xpress SARS-CoV-2' and res2.ANALYTE='SARS-CoV-2')
     )AS geneXpertSARScov2,
     
     --WGS
     (SELECT 
      wgsCovid.FINAL
      FROM RESULTS wgsCovid
      where wgsCovid.ORDNO = ot.ORDNO and wgsCovid.TESTNO = 'COVID-19 - Surveillance WGS' and wgsCovid.ANALYTE = 'Envoi' and wgsCovid.REPORTABLE = 'Y'
      AND wgsCovid.ORIGREC = ( select max(origrec) from results res2 where res2.ordno=wgsCovid.ordno and res2.testno=wgsCovid.testno AND res2.TESTNO='COVID-19 - Surveillance WGS' and res2.ANALYTE = 'Envoi')
     ) AS WGS_COVID,
     
     ( select 
      case 
        when coroLSPQ.Final > 100
        then '0'
        else  TO_CHAR(ROUND(coroLSPQ.Final))
      end
    
       from RESULTS coroLSPQ
        where coroLSPQ.ORDNO = ot.ORDNO  AND coroLSPQ.TESTNO='2019-nCoV'  AND coroLSPQ.ANALYTE='CoVN (LSPQ) Ct'  
        and coroLSPQ.ORIGREC = ( select max(origrec) from results res2 where res2.ordno=coroLSPQ.ordno and res2.testno=coroLSPQ.testno AND res2.TESTNO='2019-nCoV' and res2.ANALYTE='CoVN (LSPQ) Ct' )
     ) as CT,
     
     (
     select REPLACE(REPLACE(max(concat(concat(concat(concat(concat(adress,', '),city),', '),state),', Canada')),CHR(10),' '),CHR(13),' ') FROM RASCLIENTS where compname = rc.COMPNAME
     group by compname)as CH_ADRESS
      
 FROM 
 CENTRALRECEIVING cr 
		inner join FOLDERS f on f.FOLDERNO = cr.EXTERNAL_ID 
		inner join METADATA crm on crm.ID = cr.METADATA_GUID 
		inner join RASCLIENTS rc on rc.RASCLIENTID = cr.RASCLIENTID 
		left join PATIENTS p on p.PID = cr.PID 
		inner join ORDERS o on o.FOLDERNO = cr.EXTERNAL_ID 
		inner join ORDTASK ot on ot.ORDNO = o.ORDNO 
        
        where ot.TESTCODE IN (2685,2689,2704,2705) and cr.DATE_RECEIVED >= to_date('2020-01-01', 'YYYY-MM-DD') 
        and ( cr.PANEL_LIST like '2019-nCoV%' OR cr.PANEL_LIST like '%Soustraitant%' or cr.PANEL_LIST like '%2019-nCoV - LNM' or cr.PANEL_LIST ='GeneXpert SARS-CoV-2' 
        or cr.PANEL_LIST like '%COVID-19 - Surveillance WGS%' or cr.PANEL_LIST like '%COVID-19 - Confirmation%')
	and rc.RASCLIENTID not in ('LSPQCEC','LSPQCIC','LSPQF','LSPQP','LSPQV','LSPQ') 
        
        --AND  f.FOLDERNO = 'L00266938' AND ROWNUM < 5 order by f.FOLDERNO 
        order by f.FOLDERNO 
);








