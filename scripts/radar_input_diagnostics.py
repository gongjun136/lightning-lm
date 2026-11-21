"""逐路原始输入与组帧证据展示；不更改故障判定，不把相关性写成 DDS 根因。"""
import math


def finite(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def ids(value):
    if value in (None, 'none', ''):
        return []
    result = []
    for token in str(value).split('|'):
        n = finite(token)
        if n is None or n != int(n) or not 0 <= n <= 255:
            return []
        result.append(int(n))
    return list(dict.fromkeys(result))


def ms(value):
    n = finite(value)
    return '未记录' if n is None else f'{n * 1000:.2f} ms'


def describe_radar_inputs(evidence):
    e = evidence or {}
    if e.get('radar_rx_schema') is None:
        if e.get('radar_rx_link_status') in ('ambiguous', 'no_match'):
            reason = '符合条件的候选不唯一' if e['radar_rx_link_status'] == 'ambiguous' else '状态或帧年龄不吻合'
            return [f'逐路雷达：快照缺证（{reason}，候选 {e.get("radar_rx_link_candidates", "未知")} 条），无法归因。']
        return ['逐路雷达：此记录未含逐路接收快照，不能确定哪路导致组帧停滞（旧日志无法补算）。']
    if finite(e.get('radar_rx_schema')) != 1:
        return ['逐路雷达：快照版本不支持，无法归因。']
    if finite(e.get('radar_rx_available')) != 1:
        return ['逐路雷达：本次输入状态快照未取得，无法归因；不沿用上一次状态。']
    configured = ids(e.get('radar_rx_ids'))
    if not configured:
        return ['逐路雷达：配置来源列表缺失或无效，无法归因。']
    linked = bool(e.get('radar_rx_link_source'))
    lines = [('逐路雷达（关联证据：来自同次运行日志，非事件时间线原记录；峰值窗口为自上次快照以来）：' if linked else
              '逐路雷达（对应上方所标时间的 Safety 记录；峰值窗口为自上次快照以来）：')]
    if linked:
        lines.append(f'快照来源：同次运行日志 Safety 记录 {e.get("radar_rx_link_time", "时间未知")}，'
                     f'日志记录时间减事件时间线记录时间={ms(e.get("radar_rx_link_offset_sec"))}（关联上限 ±50 ms，非接收延迟）；'
                     '唯一候选；策略、故障状态、恢复计数一致，帧年龄在日志舍入精度内吻合；'
                     f'证据 {e["radar_rx_link_source"]}。')
    coherent = finite(e.get('radar_rx_consistent')) == 1
    after = finite(e.get('radar_rx_after_safety_sec'))
    counts = [e.get('radar_rx_safety_accepted'), e.get('radar_rx_snapshot_accepted')]
    have_counts = all(finite(v) is not None for v in counts)
    advanced = have_counts and finite(counts[1]) > finite(counts[0])
    coherent = coherent and have_counts and finite(counts[0]) == finite(counts[1])
    progress = '：未知（证据不完整）'
    if have_counts:
        progress = '已推进' if advanced else '未推进' if counts[0] == counts[1] else '回退（证据不一致）'
        progress += f'（{counts[0]}→{counts[1]}）'
    label = '快照对齐（日志内）：晚于该日志 Safety 判定' if linked else '快照对齐：晚于 Safety 判定'
    lines.append(f'{label} {ms(after)}（独立的取证上限 50 ms）；合格组帧计数{progress}。')
    if not coherent or after is None or not 0 <= after <= .05:
        recovery = '组帧推进可能已使码 2 清除，以后续 Safety 记录为准' if advanced else '是否已清除码 2，以后续 Safety 记录为准'
        lines.append('快照未对齐：以下只作补充，不认定为触发瞬间原因；' + recovery + '。')
    if (finite(e.get('lidar_frame_recovery_count')) or 0) > 0:
        lines.append('以下描述恢复确认期间的当前组帧现场，不解释码 2 最初触发的原因。')
    primary = finite(e.get('radar_rx_primary'))
    timeout = finite(e.get('radar_rx_timeout_sec'))
    missing = ids(e.get('radar_rx_pending_missing'))
    silent, delayed, rejected, unknown = [], [], [], []
    labels = {'empty': '空点云', 'invalid_stamp': '无效时间戳', 'late_bucket': '所属时间桶已排出',
              'duplicate_bucket': '同雷达同桶重复', 'tolerance': '同桶时间跨度超容差',
              'invalid': '输入无效', 'watchdog_rejected': '未通过输入组帧监视器检查',
              'accepted': '已接纳输入（不等于合格组帧）', 'unseen': '尚无回调'}
    for lidar in configured:
        prefix = f'radar_rx_{lidar}_'
        get = lambda key: e.get(prefix + key)
        count, age = finite(get('samples')), finite(get('age_sec'))
        source, arrival = finite(get('source_gap_sec')), finite(get('arrival_gap_sec'))
        state = '最后一次输入：' + labels.get(get('assembly'), '接纳状态未记录')
        stamp_state = get('stamp_state')
        if count == 0:
            silent.append(lidar)
        elif count is None or age is None:
            unknown.append(lidar)
        elif timeout is not None and age >= timeout:
            silent.append(lidar)
            state += '，此后无回调'
        if timeout is not None and source is not None and arrival is not None and 0 < source <= timeout < arrival:
            delayed.append(lidar)
            state += '，相邻回调间隔拉长而源间隔未超时'
        if get('assembly') in ('empty', 'invalid_stamp', 'late_bucket', 'duplicate_bucket', 'tolerance', 'invalid', 'watchdog_rejected'):
            rejected.append(lidar)
        if stamp_state in ('invalid', 'duplicate', 'rollback'):
            state += '，' + {'invalid': '源时间戳无效', 'duplicate': '源时间戳重复', 'rollback': '源时间戳回退'}[stamp_state]
        lag = finite(get('lag_sec'))
        if lag is not None and lag < 0:
            state += '，接收时间差为负，需核对时钟/时间戳'
        role = '主' if lidar == primary else '辅'
        window = finite(get('window_samples'))
        peak = '不适用（窗口内无回调）' if window == 0 else ms(get('max_lag_sec'))
        if window is None:
            peak += '（窗口回调数未知）'
        lines.append(f'  雷达 {lidar}（{role}，{get("topic") or "话题未记录"}）：距回调 {ms(age)}；'
                     f'最后接收时间差 {ms(lag)} / 窗口峰值 {peak}；'
                     f'源/回调间隔 {ms(source)} / {ms(arrival)}；{state}')
    if finite(e.get('radar_rx_pending')) == 1:
        present = ids(e.get('radar_rx_pending_present'))
        lines.append(f'组帧现场：最早待组帧桶已有 {present}，缺少 {missing}；'
                     f'最低来源数 {e.get("radar_rx_min_lidars", "未记录")}（且必须含主雷达）；'
                     f'同桶容差 {ms(e.get("radar_rx_tolerance_sec"))}；'
                     f'源时间推进 {ms(e.get("radar_rx_pending_progress_sec"))} / '
                     f'乱序等待窗口 {ms(e.get("radar_rx_reorder_sec"))}（按源时间推进）。')
    else:
        lines.append('组帧现场：无待组帧桶记录；不能据此认定全部雷达正常。')
    last = e.get('radar_rx_group_result')
    if last in ('missing_primary', 'stale_primary', 'watchdog_rejected'):
        reason = {'missing_primary': '缺少主雷达', 'stale_primary': '主雷达输入已过期',
                  'watchdog_rejected': '未通过组帧监视器检查'}[last]
        lines.append(f'最近一次组帧拒绝：{reason}，距该次检查 {ms(e.get("radar_rx_group_age_sec"))}。')
    if coherent and after is not None and 0 <= after <= .05:
        facts = []
        if len(silent) == len(configured):
            facts.append('所有配置雷达均未到达或回调已超时，不能锁定单一路')
            if finite(e.get('radar_rx_pending')) == 1:
                facts.append('待组帧桶为此前留下的状态，不按桶内缺路单独归因')
        elif silent:
            facts.append(f'雷达 {silent} 未到达或回调已超时')
            if primary in silent:
                facts.append('其中主雷达是合格组帧的必需来源')
        if delayed:
            facts.append(f'雷达 {delayed} 有相邻消息延迟交付线索')
        if rejected:
            facts.append(f'雷达 {rejected} 的最近输入被拒绝（原因见逐路行）')
        if unknown:
            facts.append(f'雷达 {unknown} 接收证据不完整')
        if finite(e.get('radar_rx_pending')) == 1 and len(silent) != len(configured):
            present = ids(e.get('radar_rx_pending_present'))
            minimum = finite(e.get('radar_rx_min_lidars'))
            if primary in missing:
                facts.append(f'待组帧桶缺主雷达 {int(primary)}，该桶不合格')
                if primary not in silent:
                    facts.append('主雷达仍有回调，但未进入该桶的原因未知，因果未闭合')
            elif primary in present:
                facts.append('待组帧桶含主雷达')
            else:
                facts.append('该桶主雷达是否齐备未知')
            if minimum is not None:
                enough = len(present) >= minimum
                facts.append(f'该桶来源数 {len(present)}/{int(minimum)}，' + ('已满足数量门槛' if enough else '未满足数量门槛'))
                if minimum == 1:
                    facts.append('本门槛下辅助雷达缺失单独不构成码 2')
                elif not enough:
                    facts.append('缺失来源影响该桶组帧；单个桶仍不能证明整段码 2 的唯一原因')
            else:
                facts.append('最低来源门槛未记录，不判断缺辅助雷达的影响')
        lines.append('逐路结论：' + '；'.join(facts or ['未发现足以锁定单路原因的接收证据']) + '。')
    else:
        lines.append('逐路结论：不归因（快照未对齐）。')
    lines.append('接收时间差=回调 ROS 时刻−消息时间戳；时钟同步未验证，含采样/驱动/传输/调度影响，不能直接当作 DDS 延迟。缺辅助雷达不必然触发码 2，仍取决于主雷达及实际组帧门槛。')
    return lines
