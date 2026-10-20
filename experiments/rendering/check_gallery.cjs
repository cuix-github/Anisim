// Isolated headless browser QA of the local catalog; no Maya or render process.
const fs = require('node:fs');
const path = require('node:path');
const {pathToFileURL} = require('node:url');
const assert = require('node:assert/strict');
const {chromium} = require(process.env.ANISIM_PLAYWRIGHT || 'playwright');
const root = path.resolve(__dirname, '../..');
const out = path.join(root, 'out/repo_organization_20261019');
(async () => {
  const browser = await chromium.launch({headless:true, executablePath:process.env.ANISIM_BROWSER_EXECUTABLE});
  try {
    const page = await browser.newPage({viewport:{width:1440,height:1000}});
    const errors = [], remote = [];
    page.on('pageerror',e=>errors.push(String(e)));
    page.on('request',r=>{if(/^https?:/.test(r.url()))remote.push(r.url())});
    await page.goto(pathToFileURL(path.join(root,'out/gallery/index.html')).href);
    await page.waitForFunction(()=>[...document.querySelectorAll('.preview img')].every(i=>i.complete&&i.naturalWidth>0));
    assert.equal(await page.locator('.card').count(),10);
    await page.screenshot({path:path.join(out,'gallery_desktop.png'),fullPage:true});
    await page.getByRole('button',{name:'查看 mix_u60_c40',exact:true}).click();
    assert.equal(await page.locator('dialog[open]').count(),1);
    await page.waitForFunction(()=>document.querySelector('#viewer img').naturalWidth===960);
    assert((await page.locator('#detailtitle').textContent()).includes('22 m'));
    assert.equal(await page.getByRole('link',{name:'泡沫 ASS',exact:true}).count(),1);
    await page.keyboard.press('Escape');
    await page.selectOption('#ratio','U60 / C40');
    assert.equal(await page.locator('.card').count(),1);
    await page.selectOption('#time','24.9');
    assert.equal(await page.locator('.card').count(),1);
    await page.fill('#search','nonexistent-result');
    assert.equal(await page.locator('.card').count(),0);
    await page.getByRole('button',{name:'重置',exact:true}).click();
    await page.getByRole('button',{name:'视频',exact:true}).click();
    const videos = await page.locator('.card').count(); assert(videos>0);
    await page.locator('.preview').first().click();
    await page.waitForFunction(()=>document.querySelector('video').readyState>=1,{},{timeout:15000});
    await page.getByRole('button',{name:'关闭详情'}).click();
    await page.getByRole('button',{name:'帧序列',exact:true}).click();
    const sequences=await page.locator('.card').count();assert(sequences>0);
    await page.locator('.preview').first().click();
    const options=await page.locator('.sequence option').count();assert(options>=8);
    await page.selectOption('.sequence',{index:options-1});
    await page.waitForFunction(()=>document.querySelector('#viewer img').complete&&document.querySelector('#viewer img').naturalWidth>0);
    await page.getByRole('button',{name:'关闭详情'}).click();
    await page.getByRole('button',{name:'全部档案',exact:true}).click();
    await page.selectOption('#status','unreadable');
    assert.equal(await page.locator('.card').count(),1);
    await page.getByRole('button',{name:'重置',exact:true}).click();
    await page.setViewportSize({width:390,height:844});
    await page.screenshot({path:path.join(out,'gallery_mobile.png'),fullPage:true});
    assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth));
    assert.deepEqual(errors,[]);assert.deepEqual(remote,[]);
    const result={status:'passed',selected_cards:10,video_entries:videos,sequence_cards_on_first_page:sequences,
      filter_and_search:true,detail_and_asset_links:true,video_metadata_loaded:true,frame_selection:true,
      mobile_no_horizontal_overflow:true,page_errors:errors,remote_requests:remote,simulation:false,arnold_rendering:false};
    fs.writeFileSync(path.join(out,'gallery_check.json'),JSON.stringify(result,null,2));
    console.log(JSON.stringify(result,null,2));
  } finally {await browser.close()}
})().catch(e=>{console.error(e);process.exitCode=1});
