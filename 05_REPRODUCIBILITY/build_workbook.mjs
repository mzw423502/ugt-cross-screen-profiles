import fs from 'node:fs/promises';
import path from 'node:path';
import {Workbook,SpreadsheetFile} from '@oai/artifact-tool';

// Pass the package path when running a copy of this builder with the bundled runtime.
const root=path.resolve(process.argv[2] || '.');
const data=JSON.parse(await fs.readFile(path.join(root,'05_REPRODUCIBILITY/generated/workbook_content.json'),'utf8'));
const previews=path.resolve(process.argv[3] || 'internal_checks/workbook_closeout');
await fs.mkdir(previews,{recursive:true});
const wb=Workbook.create();
for(const item of data){
  const s=wb.worksheets.add(item.name), rows=item.rows, n=rows.length, m=rows[0].length;
  s.showGridLines=false;
  const range=s.getRangeByIndexes(0,0,n,m);
  range.values=rows;
  range.format.font={name:'Arial',size:10,color:'#202020'};
  range.format.verticalAlignment='top';
  range.format.rowHeight=22;
  range.format.wrapText=true;
  for(let j=0;j<m;j++){
    const max=Math.max(...rows.map(r=>String(r[j]??'').length));
    const w=Math.min(70,Math.max(18,Math.ceil(max*1.05+2)));
    s.getRangeByIndexes(0,j,n,1).format.columnWidth=w;
    const num=rows.slice(1).map(r=>r[j]).filter(v=>v!==null);
    const numeric=num.filter(v=>typeof v==='number');
    if(numeric.length)
      s.getRangeByIndexes(1,j,n-1,1).setNumberFormat(numeric.every(Number.isInteger)?'0':'0.000');
    if((item.name==='S1_Enzyme_matches' && [2,3].includes(j)) || (item.name==='S3_Reaction_states' && j===2))
      s.getRangeByIndexes(1,j,n-1,1).setNumberFormat(item.name==='S1_Enzyme_matches'?'0.00%':'0.0%');
  }
  const h=s.getRangeByIndexes(0,0,1,m);
  h.format.fill='#e8eef2';h.format.font={name:'Arial',size:10,bold:true,color:'#152d3b'};
  h.format.borders={bottom:{style:'thin',color:'#526778'}};
  range.format.autofitRows();
  s.freezePanes.freezeRows(1);
}
wb.recalculate();
console.log((await wb.inspect({kind:'table',range:'S3_Reaction_states!A1:C5',include:'values',tableMaxRows:5,tableMaxCols:3,maxChars:1800})).ndjson);
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:10},summary:'error scan'})).ndjson);
for(const item of data){
  const rows=Math.min(item.rows.length,item.name==='S14_Biochemical_context'?8:6);
  const cols=Math.min(item.rows[0].length,8);
  const end=String.fromCharCode(64+cols)+rows;
  const blob=await wb.render({sheetName:item.name,range:'A1:'+end,scale:1.3,format:'png'});
  await fs.writeFile(path.join(previews,item.name+'.png'),new Uint8Array(await blob.arrayBuffer()));
}
const output=await SpreadsheetFile.exportXlsx(wb);
await output.save(path.join(root,'03_SUPPLEMENT/Supplementary_Data.xlsx'));
console.log('Exported Supplementary_Data.xlsx with '+data.length+' sheets');
