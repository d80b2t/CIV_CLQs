SELECT
  q.ra,q.dec,
  c.bestobjid,
  c.specobjid,
  c.ibest-q.extinction_i as i,
  c.ibesterr,
  q.extinction_i,
  q.z as zSDSS,
  c.zemHW,
  c.zemHWalt,
  c.delgi,
  so.sn_2,
  c.firstpeak,
  isnull(f.integr,-1) as integr,
  isnull(r.cps,-9) as xray,
  sl1.lineid as CIVlineID,
  sl1.restWave as CIVrestwave,
  sl1.wave as CIVwave,
  sl1.waveErr as CIVwaveErr,
  sl1.sigma as CIVsigma,
  sl1.sigmaErr as CIVsigmaErr,
  sl1.height as CIVheight,
  sl1.ew as CIVew,
  sl1.ewerr as CIVewErr,
  sl1.continuum as CIVcont,
  sl2.lineid as MgIIlineid,
  sl2.restWave as MgIIrestwave,
  sl2.wave as MgIIwave,
  sl2.waveErr as MgIIwaveErr,
  sl2.sigma as MgIIsigma,
  sl2.sigmaErr as MgIIsigmaErr,
  sl2.height as MgIIheight,
  sl2.ew as MgIIew,
sl2.ewerr as MgIIewErr,
sl2.continuum as MgIIcont, (sl1.wave/1549.06)-1 as c4z, (sl2.wave/2798.75)-1 as mg2z, ((c.zemHW)-((sl1.wave/1549.06)-1))*2.9979e5/
(1.0+c.zemHW) as c4b, ((c.zemHW)-((sl3.wave/1908.73)-1))*2.9979e5/
   (1.0+c.zemHW) as c3b
  INTO myDB.dr7qsos_c4b_HW10pub2n_big
  FROM public.gtr.DR7qsosHW10pub2  as c
  left outer join SpecPhotoAll as q on
q.specObjID = c.specObjID
left outer join SpecLine sl1 on q.specObjID =
sl1.specObjID
left outer join SpecLine sl2 on q.specObjID =
sl2.specObjID
left outer join SpecLine sl3 on q.specObjID =
sl3.specObjID
left outer join First as f on q.objID =
f.objID
left outer join Rosat as r on q.objID =
   r.objID
  left outer join SpecObjAll as so on
q.specObjID = so.specObjID
WHERE (
q.z>=1.54 AND q.z<=4.5 AND sl1.lineID=1549 AND sl2.lineID=2799 )
ORDER by q.ra
